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
ai  = 1;
cptraces = [];
for s = 1:length(subs)
    
    sub         = subs{s};  disp(['Subject ' num2str(s) ': ' sub]);
    eventid     = find(strcmp(subevents, sub));
    subE        = E(eventid); 

    M           = spm_eeg_load([Fmeeg fs sub]);
    Fs          = fsample(M);
    thesechans  = chanlabels(M);
    CPchanlab   = {'E62'}; % {'E54', 'E55', 'E80', 'E62'};
    C3chanlab   = {'E37'};
    C4chanlab   = {'E105'};

    % Find IDs for CP channels
    %----------------------------------------------------------------------
    i = 1;
    for c = 1:length(CPchanlab)
        match   = find(strcmp(thesechans, CPchanlab{c}));
        if ~isempty(match)
            CPid(i)     = match;
            i           = i + 1;
        end
    end
    
    % Find IDs for C3 channels
    %----------------------------------------------------------------------
    i = 1;
	for c = 1:length(C3chanlab)
        match   = find(strcmp(thesechans, C3chanlab{c}));
        if ~isempty(match)
            C3id(i)     = match;
            i           = i + 1;
        end
    end
    
    % Find IDs for C4 channels
    %----------------------------------------------------------------------
	for c = 1:length(CPchanlab)
        match   = find(strcmp(thesechans, CPchanlab{c}));
        if ~isempty(match)
            CPid(i)     = match;
            i           = i + 1;
        end
    end
    
    % Collate time traces for CP channel
    %======================================================================
    for m = 1:size(M,3)
        cptraces(ai,:) = mean(M(CPid,:,m),1);
        ai  = ai + 1;
    end
    
    % Sliding window estimator of beta frequency for C3/4 chanels
    %======================================================================
    for m = 1:size(M,3)
        win     = 0.4  * Fs;    % length of window in seconds
        step    = 0.05 * Fs;    % step size in seconds
        lseg    = nsamples(M);  
        
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
            dat     = M(C3id,thiswin,:);  % channels * time * trials
            fdat   = fft(dat, [], 2);
            fbystep(i,:,:) = mean(fdat(:,Nhalf,:),1);
        end
        
        
        
    end
end

%%
rctm = [E.rctm];

Fs          = fsample(M);
locked      = [];

for c = 1:size(cptraces,1)
    smpl = fix((rctm(c) + 500) / 1000 * Fs);
    start   = fix(smpl - .7*Fs);
    stop    = fix(smpl + .1*Fs)-1;
    
    trace   = cptraces(c,start:stop);
    if sum(abs(trace) > 100) > 1, trace = zeros(1,length(trace)); end
    trace   = trace - mean(trace(1:125));
    locked(c,:)     = trace;
    
end
plot(mean(locked,1)); hold on
