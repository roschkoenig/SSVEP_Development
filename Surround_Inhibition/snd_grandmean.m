%% Housekeeping
%==========================================================================
clear all
D   = snd_housekeeping;
fs  = filesep;

Fscripts    = D.Fscripts;
Fdata       = D.Fdata;
Fgenex      = D.Fgenex;
Fanalysis   = D.Fanalysis;
Fmeeg       = D.Fmeeg;
clear D 

% Load list of subjects with MEEG files available
%--------------------------------------------------------------------------
schar       = spm_select('List', [Fanalysis fs 'MEEG'], '^A.*.mat');
dotpos      = find(schar(1,:) == '.');  Apos = find(schar(1,:) == 'A');
subs        = cellstr(schar(:,Apos:dotpos-1));
%%
for s = 1:length(subs)
    disp(['Subject No: ' num2str(s)]); 
    sub         = subs{s};
    
    job{1}.spm.meeg.averaging.average.D = {[Fmeeg fs sub '.mat']};
    job{1}.spm.meeg.averaging.average.userobust.standard = false;
    job{1}.spm.meeg.averaging.average.plv = false;
    job{1}.spm.meeg.averaging.average.prefix = 'm';
    spm_jobman('run', job);
    clear job
    
end
%%
meanfiles = cellstr(spm_select('FPList', Fmeeg, '^m.*.mat'));
job{1}.spm.meeg.averaging.grandmean.D = meanfiles;
job{1}.spm.meeg.averaging.grandmean.outfile = 'grandmean';
job{1}.spm.meeg.averaging.grandmean.weighted = 1;
spm_jobman('run', job);
clear job
