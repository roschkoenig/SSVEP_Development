%% Housekeeping
%==========================================================================
clear all
D   = snd_housekeeping;
fs  = filesep;
Fbase       = D.Fbase;
Fscripts    = D.Fscripts;
Fdata       = D.Fdata;
Fgenex      = D.Fgenex;
Fanalysis   = D.Fanalysis;
Fmeeg       = D.Fmeeg;
clear D 

extract_trigger     = 0;
cond_id             = 1;

% Define subjects for subject loop
%--------------------------------------------------------------------------
subs    = cellstr(spm_select('List', [Fdata fs 'EEG'], 'dir', '^A'));

% Extract trigger data from EEG files
%==========================================================================
if extract_trigger
for sb = 1:length(subs)
    clearvars -except Fscripts Fdata Fgenex Fanalysis Fmeeg subs sb fs alles cond_id
    
    disp(['Subject: ' num2str(sb)]);
    disp('Identifying datafiles');
    
    sub     = subs{sb};
    Feeg	= [Fdata fs 'EEG' fs sub fs 'EEG' fs 'preprocessed' fs 'mat_format'];
    datafiles = cellstr(spm_select('FPList', Feeg, '^(?!T).*\.mat'));
    
    clear trg
    for d = 1:length(datafiles)
        try 
        data    = load(datafiles{d}); 
        if ~isempty(data.EEG.event)
            trg(d)      = str2double(data.EEG.event(1).type);
        end
        end
        clear data;
    end
    save([Feeg fs 'Triggers'], 'trg');
    clear d
end


%% Load relevant datasets (triggers 94, 95, 96)
%--------------------------------------------------------------------------
for sb = 1:length(subs)
	disp(['Subject: ' num2str(sb)]);
    disp('Loading trigger timings');
    
	sub         = subs{sb};
    Feeg        = [Fdata fs 'EEG' fs sub fs 'EEG' fs 'preprocessed' fs 'mat_format'];
    datafiles   = cellstr(spm_select('FPList', Feeg, '^(?!T).*\.mat'));
    load([Feeg fs 'Triggers']);
    
    c = 0;
    fileid  = [];
    for t = [94 95 96]
        tid = find(trg == t);
        if ~isempty(tid)
        for tt = tid
            c               = c + 1;
            fileid          = [fileid, tt];
        end
        end
    end
    fileids{sb} = fileid;
    save([Fanalysis fs 'FileIDs'], 'fileids');
end
end
    
%% Identify conditions
%==========================================================================
% The continuous task will be broken up into baseline and target, left and
% right, as well as correct and incorrect conditions (i.e. 8 conditions). 

if cond_id

load([Fanalysis fs 'FileIDs']); 
allcount = 0;

for sb = 1:length(subs)
	disp(['Subject: ' num2str(sb)]);
    disp('Epoching by condition');
    
	sub         = subs{sb};
    Feeg        = [Fdata fs 'EEG' fs sub fs 'EEG' fs 'preprocessed' fs 'mat_format'];
    datafiles   = cellstr(spm_select('FPList', Feeg, '^(?!T).*\.mat'));
    load([Feeg fs 'Triggers']);

    % Loop through all relevant epochs
    %----------------------------------------------------------------------
    if ~isempty(fileids{sb})
    for d = 1:length(fileids{sb})
        
        clear D EEG
        events = []; latency = [];
        load(datafiles{fileids{sb}(d)}, 'EEG'); 
        D = EEG;
        
        for o = 1:length(D.event)
            events(o)   = str2double(D.event(o).type);
            latency(o)  = D.event(o).latency;  
        end

        % Find event markers for the start event (marker == 5);
        %------------------------------------------------------------------
        onsets = find(events == 5);
        i  = 0;
        
        % Loop through all trials to identify decisions and reaction times
        %------------------------------------------------------------------
        clear decision decevent trialtype trialevent
        corr_ons = [];
        
        for o = 1:length(onsets)-1
        if length(onsets(o) : onsets(o+1)) == 4                             % only include correct no of triggers
            i               = i+1;
            corr_ons(i)     = onsets(o);

            for s = onsets(o):onsets(o+1)
                switch events(s) 
                    case 12                                                 % left button
                        decision(i) = 'L';
                        decevent(i) = s;  
                    case 13                                                 % right button
                        decision(i) = 'R';
                        decevent(i) = s;
                    case 8 
                        trialtype(i) = 'L';                             % left target
                        trialevent(i) = s;
                    case 9 
                        trialtype(i) = 'R';                             % right target
                        trialevent(i) = s;
                end
            end
            
        end     % If-conditional to only identify sensible trials
        end     % For-loop through each onset marker 
        
        for k = 1:length(corr_ons)
        if  trialevent(k) ~= 0 && decevent(k) ~= 0 && [latency(decevent(k)) - latency(trialevent(k)) > 400]
            
            allcount = allcount + 1;
            
            E(allcount).sub     = sub;
            E(allcount).fname   = datafiles{fileids{sb}(d)};
            
            E(allcount).side    = trialtype(k);
            E(allcount).resp    = decision(k);
            E(allcount).corr    = trialtype(k) == decision(k); 
            E(allcount).rctm    = latency(decevent(k)) - latency(trialevent(k)); 
            
            E(allcount).labl    = [E(allcount).side E(allcount).resp];
            E(allcount).run     = d;
            
            % Convert timing details into sample indices
            %--------------------------------------------------------------
            eventtime           = latency(trialevent(k));
            E(allcount).time    = eventtime;
            E(allcount).onsmpl  = nearest(D.times, eventtime);
            
        end
        end
        
    end     % For-loop through individual files
    end     % If-conditional on relevant recordings identified
end  

save([Fanalysis fs 'Events'], 'E');

end