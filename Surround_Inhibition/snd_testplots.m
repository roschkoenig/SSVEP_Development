%% Housekeeping
%==========================================================================
D   = snd_housekeeping;
fs  = filesep;

Fscripts    = D.Fscripts;
Fdata       = D.Fdata;
Fgenex      = D.Fgenex;
Fanalaysis  = D.Fanalysis;
Fmeeg       = D.Fmeeg;
clear D 

subs        = {'A00051826', 'A00054894'};
sub         = subs{2};
Feeg        = [Fdata fs 'EEG' fs sub fs 'EEG' fs 'preprocessed' fs 'mat_format'];

% Load MEEG object and extract trial info
%--------------------------------------------------------------------------
M           = spm_eeg_load([Fmeeg fs sub]);
Fs          = fsample(M);
cnds        = condlist(M);
trialcond   = conditions(M);
lbls        = chanlabels(M);


% Define conditions of interest
%--------------------------------------------------------------------------
ci{1}      = find(strcmp(cnds, 'RR'));   % BL for CR response to RT stim
ci{2}      = find(strcmp(cnds, 'LL'));

for cc = 1:length(ci)
    
cond_id = find(ismember(trialcond, cnds(ci{cc})));
chnl    = find(strcmp(lbls, 'E75'));

for c = 1:length(cond_id)
    dat{cc}(c,:)            = M(chnl,:,cond_id(c));
    [wt{cc}(c,:,:),fax]     = cwt(dat{cc}(c,:), Fs); 
    tft                     = fft(dat{cc}(c,:));
    ft{cc}(c,:)             = abs(tft(1:fix(length(tft)/2)));
    clear tft
end

subplot(2,length(ci),cc)
[val loc] = min(abs(fax - 20));
imagesc(squeeze(abs(mean(wt{cc},1)))); hold on

[val loc] = min(abs(fax - 20));
plot([0 length(dat{cc}(1,:))], [loc loc], 'Linewidth', 1, 'color', 'w');

[val loc] = min(abs(fax - 25));
plot([0 length(dat{cc}(1,:))], [loc loc], 'Linewidth', 1, 'color', 'r');

subplot(2, length(ci), cc + length(ci))
plot(squeeze(abs(mean(ft{cc},1))));
end
%%
frres = fix(length(dat{1})/2);
maxfr = fix(ceil(Fs)/2);
minfr = 1/length(dat{1}) / ceil(Fs);
frax  = linspace(minfr, maxfr, frres);

plot(frax, log(mean(ft{1},1))); hold on
plot(frax, log(mean(ft{2},1)))
xlim([1 50]);
