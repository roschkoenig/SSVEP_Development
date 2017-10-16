%% Housekeeping
%==========================================================================
clear all
D   = snd_housekeeping('surround');
fs  = filesep;

Fscripts    = D.Fscripts;
Fdata       = D.Fdata;
Fgenex      = D.Fgenex;
Fanalysis   = D.Fanalysis;
Fmeeg       = D.Fmeeg;
clear D

% Load list of subjects with MEEG files
%--------------------------------------------------------------------------
load([Fanalysis fs 'SPM' fs 'Subjectlist']);
subs    = slist;


%% Loop through individual subjects
%==========================================================================
count       = 0;
for s = 1:length(subs)
sub         = subs{s};  disp(['Subject ' num2str(s) ': ' sub]);
Feeg        = [Fdata fs 'EEG' fs sub fs 'EEG' fs 'preprocessed' fs 'mat_format'];

M           = spm_eeg_load([Fmeeg fs 'ra' sub]);
Fs          = fsample(M); 

cnds        = condlist(M);
trialcond   = conditions(M);
lbls        = chanlabels(M);

% Define groups of channels
%==========================================================================
chanlabs    = chanlabels(M);

% Oz electrode group
%--------------------------------------------------------------------------
Ozchans     = {'E75', 'E71', 'E70', 'E74', 'E81', 'E82', 'E83', 'E76'};
Ozid        = [];
i           = 0;

for o = 1:length(Ozchans)
    thischan    = Ozchans(o);
    chanid      = find(strcmp(chanlabs, thischan));
    if ~isempty(chanid)
    if [isempty(badchannels(M)) | isempty(find(badchannels(M) == chanid))] 
        i       = i + 1;
        Ozid(i) = chanid;
    end
    end
end

% Estimating SSVEP from Oz electrode group
%==========================================================================
% Sliding window fourier transform
%--------------------------------------------------------------------------
win     = 0.4  * Fs;    % length of window in seconds
step    = 0.05 * Fs;    % step size in seconds
lseg    = nsamples(M);  
clear fdat

% Prepare the storage variable for fourier transforms (fbystep)
%--------------------------------------------------------------------------
% i           = 0;
% lfinal      = length(1:step:(lseg-win));
% [pow fax]   = sri_fourier(rand(win,1), Fs);

% Prepare the timeaxis associated with sliding window (timstpsec)
%--------------------------------------------------------------------------
timesec     = time(M);
% timstpsec   = linspace(timesec(1) + 0.5*win/Fs, timesec(end) - 0.5*win/Fs, lfinal);

% % Sliding window estimator
% %--------------------------------------------------------------------------
% for ww = 1 : step : (lseg-win)
%     i = i + 1;
%     thiswin         = ww: ww+win - 1;
%     dat             = M(Ozid,thiswin,:);                                    % channels * time * trials
%     clear allfdat
%     for di = 1:size(dat,3)
%         d = squeeze(dat(:,:,di));
%         [fdat]          = sri_fourier(d', Fs);                              % freq * channels
%         allfdat(:,:,di) = fdat';                                            % channels * freq * trials
%     end
%     fbystep(i,:,:)  = mean(allfdat,1);                                      % average over channels
% end

% Extract SSVEP measures
%==========================================================================
% For each (selection) event, select the appropriate SSVEP freq
clear ssvep conds ssvep_norm
conds   = conditions(M);
for e = 1:size(M,3)   
    bldat           = M(Ozid, find(timesec < 0), :); 
%     [blpow blfreq]  = sri_fourier(
    
%     % Get baseline average 
%     %----------------------------------------------------------------------
%     blid    = find(timstpsec < 0);
%     blssvep = mean(abs(fbystep(blid,fid,e)),1);
    
    bgpower         = mean(abs(fbystep(:,[fid - 1, fid + 2],e)), 2);
    ssvep_norm(:,e) = (abs(fbystep(:,fid,e))) - bgpower;
    ssvep(:,e)      = (abs(fbystep(:,fid,e)));
end

% Calculate SSVEP over the whole time window
%--------------------------------------------------------------------------
try
dat     = M(Ozid, 251:750, :);  % channels * time * trials
clear chanvep
for d = 1:size(dat,1)
    [pow freq] = sri_fourier(squeeze(dat(d,:,:)), Fs);
    chanvep(d,:,:)  = pow;
end
long.spectrum   = squeeze(mean(chanvep,1));
long.freq       = freq;

count           = count + 1;
S(s).name       = sub;
S(s).long       = long;
S(s).ssvep      = ssvep;
S(s).ssvep_norm = ssvep_norm;
S(s).fourier    = abs(fbystep);
S(s).fax        = fax;
S(s).timesteps  = timstpsec;
S(s).conds      = conds;
end

end

%% Plot Spectrum by age
%--------------------------------------------------------------------------

% Extract age information
%--------------------------------------------------------------------------
pheno   = csv2cell([Fdata fs 'MIPDB' fs 'MIPDB_PublicFile.csv'], 'fromfile');

[r agec] = find(~cellfun(@isempty, strfind(pheno, 'Age')));
agecells = pheno(2:end, agec);
for a = 1:length(agecells)
    ages(a) = str2double(agecells{a});
end

[r subc] = find(~cellfun(@isempty, strfind(pheno, 'ID')));
subcells = pheno(2:end, subc);

[r subc] = find(~cellfun(@isempty, strfind(pheno, 'DX_Status')));
diagnosis = pheno(2:end, subc);


% Extract individual trials and corresponding ages
%--------------------------------------------------------------------------
clear shortspec agelist ssvep
count = 0;
for s = 1:length(S)
    sub     = S(s).name;
    subi    = find(strcmp(subcells, sub));
    
    age(s)  = ages(subi);
    dx(s)   = diagnosis{subi};
    
    for si = 1:size(S(s).long.spectrum, 2)
        if strcmp(S(s).conds{si}, 'F3B0')
            fax  = S(s).long.freq;
            Hz25 = nearest(fax, 25);
            faxi = intersect(find(S(s).long.freq > 1), find(S(s).long.freq < 40));
            
            count = count + 1;
            shortspec(count,:)  = S(s).long.spectrum(faxi, si);
            agelist(count)      = age(s);
            ssvep(count)        = S(s).long.spectrum(Hz25, si) - S(s).long.spectrum(Hz25-1, si);
        end
    end
end

% Sort extracted spectra
%--------------------------------------------------------------------------
[sorted sorting] = sort(agelist);
sortedspec          = shortspec(sorting,:);

win     = 40;
step    = 1;

count = 0; 
clear smoothspec
for i = 1:1:(size(sortedspec,1)-win)
    count = count + 1;
    smoothspec(count,:) = mean(sortedspec(i:i+win, :), 1);
end


%%
cols = flip(cbrewer('div', 'Spectral', size(smoothspec,1)));
for s = 1:size(smoothspec,1)
    plot(fax(faxi), smoothspec(s,:), 'color', cols(s,:)); hold on
end

%%




for i = 0:3

spectra = [];
for s = 1:length(S)
nobg = find(~cellfun(@isempty, strfind(S(s).conds, ['F' num2str(i)])));
spectra = [spectra, S(s).long.spectrum(:,nobg)];
end

plot(S(1).long.freq, log(mean(spectra,2))); hold on
xlim([1 60]); 
end

legend({'1', '2', '3', '4'});
%%
count = 0;
clear allvep allfg allbg thisvep Fi Bi
C = cell(4,2);
for s = 1:length(S)
    thistime    = find(S(s).timesteps > 0);
    thisvep     = mean(S(s).ssvep_norm(thistime, :),1);
    for t = 1:length(thisvep)
        count           = count + 1;
        allvep(count)   = thisvep(t);
        Fi(count)       = str2double(S(s).conds{t}(2)) + 1;
        Bi(count)       = str2double(S(s).conds{t}(4)) + 1;
        
        C{Fi(count),Bi(count)} = [C{Fi(count),Bi(count)}, allvep(count)];
    end
end

for f = 1:size(C,1)
for b = 1:size(C,2)
    SSVEP(f,b) = mean(C{f,b});
end
end

plot(SSVEP)
legend('No background', 'Background');


