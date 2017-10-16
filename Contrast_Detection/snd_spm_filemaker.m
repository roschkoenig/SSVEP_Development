%% Housekeeping
%==========================================================================
clear all
D   = snd_housekeeping('contrast');
fs  = filesep;
Fbase       = D.Fbase;
Fscripts    = D.Fscripts;
Fdata       = D.Fdata;
Fgenex      = D.Fgenex;
Fanalysis   = D.Fanalysis;
Fmeeg       = D.Fmeeg;

EEGpath     = 'EEG old';
Feeg        = [Fdata fs EEGpath];

% Define subjects for subject loop
%--------------------------------------------------------------------------
subs    = cellstr(spm_select('List', [Fdata fs EEGpath], 'dir', '^A'));

% Loop through subjects 
%==========================================================================
for sb = 1:length(subs)
    disp(['Subject: ' num2str(sb)]);
    sub         = subs{sb};
    Feeg        = [Feeg fs sub fs 'EEG' fs 'preprocessed' fs 'mat_format'];
    
    if ~isempty(eventid)   
        
        % Make MEEG File
        %==================================================================
        clear ftdata allnames
        
        subE    = E(eventid); 
        for s = 1:length(subE)
            allnames{s} = subE(s).fname;
        end

        names   = unique(allnames, 'stable');
        c       = 1;        
        
        for d = 1:length(names)            
            fname = names{d};
            load(fname, 'EEG');     D = EEG;    clear EEG; 
            
            Fs      = D.srate;                                              % Sampling frequency
            prew    = 0.5;                                                  % prestimulus window in seconds
            win     = 2.5;                                                  % total length of window of interest in seconds
        
            eid = find(strcmp(allnames, fname));
            for ee = 1:length(eid)
                clear Ee
                e   = eid(ee);
                Ee  = subE(e);
                
                start = Ee.onsmpl - (prew * Fs);
                stop  = Ee.onsmpl + win * Fs - 1;
                
                ftdata.trial{c} = D.data(:, start:stop);
                ftdata.time{c}  = linspace(-prew, win - 1/Fs, length(ftdata.trial{c})); 
                c               = c+1;
            end
            ftdata.label    = {D.chanlocs.labels};
            fname           = [Fmeeg fs sub '.mat'];
        end
        
        % Rereference to average reference
        %------------------------------------------------------------------
        for t = 1:length(ftdata.trial)
            ftdata.trial{t} = ft_preproc_rereference(ftdata.trial{t});
%             ftdata.trial{t} = ft_preproc_bandpassfilter(ftdata.trial{t}, Fs, [1 80]);
        end
        
        % Generate MEEG object
        %--------------------------------------------------------------------------
        M = spm_eeg_ft2spm(ftdata, fname);
        M = type(M, 'single'); 

        for c = 1:length(subE)
            M = conditions(M, c, subE(c).labl);
        end
        save(M)

        % Load EEG sensor locations and project to scalp maps
        %------------------------------------------------------------------
        S           = [];
        S.task      = 'loadeegsens';
        S.source    = 'locfile';
        S.sensfile  = [Fdata fs 'EEG' fs 'GSN_HydroCel_129.sfp'];
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

%         job{1}.spm.meeg.preproc.artefact.D = {fname};
%         job{1}.spm.meeg.preproc.artefact.mode = 'reject';
%         job{1}.spm.meeg.preproc.artefact.badchanthresh = 0.2;
%         job{1}.spm.meeg.preproc.artefact.append = true;
%         job{1}.spm.meeg.preproc.artefact.methods.channels{1}.all = 'all';
%         job{1}.spm.meeg.preproc.artefact.methods.fun.peak2peak.threshold = 150;
%         job{1}.spm.meeg.preproc.artefact.prefix = 'a';
%         
%         job{2}.spm.meeg.preproc.remove.D = {[Fmeeg fs 'a' sub '.mat'];};
%         job{2}.spm.meeg.preproc.remove.prefix = 'r';
%         
%         spm_jobman('run', job);
%         clear job

    end

end
