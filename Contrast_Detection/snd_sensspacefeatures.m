%% Housekeeping
%==========================================================================
clear all
D   = snd_housekeeping;
fs  = filesep;

Fscripts    = D.Fscripts;
Fdata       = D.Fdata;
Fgenex      = D.Fgenex;
Fanalysis  = D.Fanalysis;
Fmeeg       = D.Fmeeg;
clear D

% Load list of subjects with MEEG files
%--------------------------------------------------------------------------
subs        = spm_select('List', [Fanalysis fs 'MEEG'], '^A.*.mat');
subs        = cellstr(subs(:,1:9));

% Load events file
%--------------------------------------------------------------------------
load([Fanalysis fs 'Events']);
subevents = cellstr(vertcat(E.sub));


% Loop through individual subjects
%==========================================================================
for s = 1:length(subs)
sub         = subs{s};  disp(['Subject ' num2str(s) ': ' sub]);
Feeg        = [Fdata fs 'EEG' fs sub fs 'EEG' fs 'preprocessed' fs 'mat_format'];
eventid     = find(strcmp(subevents, sub));
subE        = E(eventid); 

M           = spm_eeg_load([Fmeeg fs sub]);
Fs          = fsample(M); 

cnds        = condlist(M);
trialcond   = conditions(M);
lbls        = chanlabels(M);

% Define groups of channels
%==========================================================================
chanlabs    = chanlabels(M);

% Oz electrode group
%--------------------------------------------------------------------------
Ozchans     = {'E75', 'E71', 'E70', 'E74', 'E82', 'E83', 'E76'};
Ozid        = [];
i           = 0;

for o = 1:length(Ozchans)
    thischan    = Ozchans(o);
    chanid      = find(strcmp(chanlabs, thischan));
    if [isempty(badchannels(M)) | ~isempty(find(badchannels(M) == chanid))] & ~isempty(chanid)
        i       = i + 1;
        Ozid(i) = chanid;
    end
end

% CPP electrode group
%--------------------------------------------------------------------------
CPchans     = {'E54', 'E55', 'E80'}; %%{'E62', 'E67', 'E77'}; %, 'E61', 'E68', 'E79', 'E8cc0'}; 
CPid        = [];
i           = 0;

for c = 1:length(CPchans)
    thischan    = CPchans(c);
    chanid      = find(strcmp(chanlabs, thischan));
    if [isempty(badchannels(M)) | ~isempty(find(badchannels(M) == chanid))] & ~isempty(chanid)
        i       = i + 1;
        CPid(i) = chanid;
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

% Define frequency axis for conversion to Hz
%--------------------------------------------------------------------------
Nhalf       = 1:ceil(win/2);
fax         = linspace(0,floor(Fs/2), Nhalf(end));

% Prepare the storage variable for fourier transforms (fbystep)
%--------------------------------------------------------------------------
i           = 0;
lfinal      = length(1:step:(lseg-win));
fbystep     = zeros(lfinal,length(Nhalf),size(M,3));    % steps * freq * trials

% Prepare the timeaxis associated with sliding window (timstpsec)
%--------------------------------------------------------------------------
timesec     = time(M);
timstpsec   = linspace(timesec(1) + 0.5*win/Fs, timesec(end) - 0.5*win/Fs, lfinal);

% Sliding window estimator
%--------------------------------------------------------------------------
for ww = 1 : step : (lseg-win)
    i = i + 1;
    thiswin = ww: ww+win - 1;
    dat     = M(Ozid,thiswin,:);  % channels * time * trials
    fdat   = fft(dat, [], 2);
    fbystep(i,:,:) = mean(fdat(:,Nhalf,:),1);
end

% Extract SSVEP measures
%==========================================================================
% For each (selection) event, select the appropriate SSVEP freq
clear ssvep
for e = 1:length(subE)
    
    % Identify which side the target is (to define frequency of interest)
    %----------------------------------------------------------------------
    side    = subE(e).side;
    clear thisfreq
    if side == 'L', thisfreq = 20; elseif side == 'R', thisfreq = 25; end
    fid     = nearest(fax, thisfreq);
    
    % Get baseline average 
    %----------------------------------------------------------------------
    blid    = find(timstpsec < 0);
    blssvep = mean(abs(fbystep(blid,fid,e)),1);
    
    ssvep(:,e) = (abs(fbystep(:,fid,e))) / blssvep;
end

% Estimating ERP from CPP electrode group
%==========================================================================
cptrace = squeeze(mean(M(CPid,:,:), 1));
cptrace = cptrace';
cptrace = ft_preproc_bandpassfilter(cptrace, Fs, [1 5], 4);

% Get baseline average 
%----------------------------------------------------------------------
blid    = find(time(M) < 0);
blcp    = mean(cptrace(blid));

cptrace = cptrace - blcp;


S(s).ssvep      = ssvep;
S(s).fourier    = abs(fbystep);
S(s).fax        = fax;
S(s).timesteps  = timstpsec;
S(s).cptrace    = cptrace;
S(s).events     = subE;

end

save([Fanalysis fs 'Fourier_summaries'], 'S');


