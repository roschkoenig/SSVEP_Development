%% Housekeeping
%==========================================================================
clear all
D   = snd_housekeeping('surround');
fs  = filesep;

Fbase       = D.Fbase;
Fscripts    = D.Fscripts;
Fdata       = D.Fdata;
Fgenex      = D.Fgenex;
Fanalysis   = D.Fanalysis;
Fmeeg       = D.Fmeeg;

EEGpath     = 'MIPDB';    % 'MIPDB' or 'HBN'
writeover   = 1;

% Define subjects for subject loop
%--------------------------------------------------------------------------
subs    = cellstr(spm_select('List', [Fdata fs EEGpath], 'dir', '^A'));

% Loop through subjects 
%==========================================================================
clear alldat 
i = 0;

for sb = 1:length(subs)

    disp(['Subject: ' num2str(sb)])
    sub         = subs{sb};
    Fpre        = [Fdata fs EEGpath fs sub fs 'EEG' fs 'preprocessed'];
    Feeg        = [Fpre fs 'mat_format'];
    
    if writeover || exist([Fmeeg fs sub '.mat'], 'file') == 0
        datafiles   = cellstr(spm_select('FPList', Feeg, '^(?!T).*\.mat'));

        clear trg
        for d = 1:length(datafiles)
            try 
                data    = load(datafiles{d});                 
                if ~isempty(data.EEG.event)
                    trg(d)      = str2double(data.EEG.event(1).type);
                else
                    trg(d) = 0; 
                end
            catch
                trg(d) = 0;
            end
            clear data;
        end
        if length(trg) == length(datafiles); disp('All Triggers found');
        else warning('Something is wrong with the file list'); end
        clear d

        % Make MEEG File
        %==================================================================
        targettrg = [93 97];    behvnames = {'SurroundSupp_Block1'; 'SurroundSupp_Block2'};
        behvfiles   = cellstr(spm_select('FPList', Fpre, '^*\.mat$'));
        allcount    = 0;
        
        if length(behvfiles) == 2
        for t = 1:length(targettrg) 

            % Identify the correct  file to load
            %------------------------------------------------------------------
            clear D
            tid = find(trg == targettrg(t));
            
            thisbehv    = find(~cellfun(@isempty, strfind(behvfiles, behvnames{t})));
            disp(['Loading ' behvnames{thisbehv}]);
            B = load(behvfiles{thisbehv});
            
            % Identify the right datafile to match the behavioural files
            %--------------------------------------------------------------
            clear eventid 
            shouldIload = 0;
            
            for tt = 1:length(tid)
                
                load(datafiles{tid(tt)}, 'EEG');    D = EEG;    clear EEG;
                for e = 1:length(D.event)
                    eventid(e) = str2double(D.event(e).type);
                end
                
                % Check we have the right number of onset times
                %----------------------------------------------------------
                if length(find(eventid == 8)) ~= length(B.TargOnT)
                    shouldIload(tt) = 0;
                else 
                    shouldIload(tt) = 1;
                end
                
                clear D
                
            end     % Looping through multiple identical triggers
            
            if sum(shouldIload) == 1
                disp('Found a matching datafile');
                load(datafiles{tid(find(shouldIload))}, 'EEG');    D = EEG;    clear EEG;

                locs    = {D.chanlocs.labels};
                Fs      = D.srate;                                              % Sampling frequency
                prew    = 0.5;                                                  % prestimulus window in seconds
                win     = 2.4;                                                  % total length of window of interest in seconds
                D.data  = ft_preproc_bandpassfilter(D.data, Fs, [.5 80]);


                % Identify all onset times and store in onsmpl (i.e. onset sample)
                %------------------------------------------------------------------
                count = 0; 
                clear E eventid
                for e = 1:length(D.event)
                    eventid(e) = str2double(D.event(e).type);
                    if eventid(e) == 8
                        count = count + 1;
                        E(count).onsmpl = D.event(e).latency;
                    end
                end

                % Extract condition labels for each onset time
                %------------------------------------------------------------------
                numfg = ceil(B.CNTcon * 3);
                for o = 1:length(E)
                    E(o).bg = ['B' num2str(B.BGcon(o))];
                    E(o).fg = ['F' num2str(numfg(o))];
                end

                % Collate data required for the MEEG object
                %------------------------------------------------------------------
                for o = 1:length(E)
                    e = E(o);
                    start = e.onsmpl - (prew * Fs);
                    stop  = e.onsmpl + win * Fs - 1;

                    allcount                = allcount + 1;
                    ftdata.trial{allcount}  = D.data(:, start:stop);
                    ftdata.time{allcount}   = linspace(-prew, win - 1/Fs, length(ftdata.trial{allcount})); 
                    cond{allcount}          = [E(o).fg E(o).bg];
                    
                end
                
            end     % Condition that only 1 datafiles matches behaviour
        end     % Loop through different triggers
        end     % Check that there are some behavioural files
        
        if allcount == 0
            disp(['Subject ' num2str(sb) ' (' sub ') does not have any data']);
            
        else
            i = i+1;
            alldat(i,:) = ftdata.trial{1}(1,:);
            
            ftdata.label    = locs;
            fname           = [Fmeeg fs sub '.mat'];           
    
            % Rereference to average reference
            %------------------------------------------------------------------
            for t = 1:length(ftdata.trial)
                ftdata.trial{t} = ft_preproc_rereference(ftdata.trial{t});
            end
    
            % Generate MEEG object
            %--------------------------------------------------------------------------
            M = spm_eeg_ft2spm(ftdata, fname);
            M = type(M, 'single'); 
    
            for c = 1:length(cond)
                M = conditions(M, c, cond{c});
            end
            save(M)
    
            % Load EEG sensor locations and project to scalp maps
            %------------------------------------------------------------------
            S           = [];
            S.task      = 'loadeegsens';
            S.source    = 'locfile';
            S.sensfile  = [Fdata fs EEGpath fs 'GSN_HydroCel_129.sfp'];
            S.D         = M;
            M           = spm_eeg_prep(S);
            save(M);
    
            S           = [];
            S.task      = 'project3D';
            S.modality  = 'EEG';
            S.D         = M;
            M           = spm_eeg_prep(S);
            save(M);
    
            % Compute leadfields in preparation for model inversion
            %======================================================================
            clear job
            job{1}.spm.meeg.source.headmodel.D = {fname};
            job{1}.spm.meeg.source.headmodel.val = 1;
            job{1}.spm.meeg.source.headmodel.comment = '';
            job{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1;
            job{1}.spm.meeg.source.headmodel.meshing.meshres = 3;
            job{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'FidNz';
            job{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
            job{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'FidT9';
            job{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'lpa';
            job{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'FidT10';
            job{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'rpa';
            job{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
            job{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
            job{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
            spm_jobman('run', job);
            clear job
    
            M = spm_eeg_inv_forward(fname);
            save(M)
    
            job{1}.spm.meeg.preproc.artefact.D = {fname};
            job{1}.spm.meeg.preproc.artefact.mode = 'reject';
            job{1}.spm.meeg.preproc.artefact.badchanthresh = 0.2;
            job{1}.spm.meeg.preproc.artefact.append = true;
            job{1}.spm.meeg.preproc.artefact.methods.channels{1}.all = 'all';
            job{1}.spm.meeg.preproc.artefact.methods.fun.peak2peak.threshold = 150;
            job{1}.spm.meeg.preproc.artefact.prefix = 'a';
    
            job{2}.spm.meeg.preproc.remove.D = {[Fmeeg fs 'a' sub '.mat'];};
            job{2}.spm.meeg.preproc.remove.prefix = 'r';
    
            spm_jobman('run', job);
            clear job
        end
    end     % if condition checking whether there are any data
end

