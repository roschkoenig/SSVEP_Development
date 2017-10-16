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
clear D 

%% Define subjects for subject loop
%--------------------------------------------------------------------------
subs    = cellstr(spm_select('List', [Fdata fs 'EEG'], 'dir', '^A'));
origpath    = '/Volumes/ANRI/DevEEG';

for s = 1:length(subs)
    sub = subs{s};
    f1 = [sub '_SurroundSupp_Block1.mat'];
    f2 = [sub '_SurroundSupp_Block2.mat'];
    
    Fpre	= [Fdata fs 'EEG' fs sub fs 'EEG' fs 'preprocessed'];
    Forig   = [origpath fs sub fs 'Behavioral' fs 'mat'];
    try
        copyfile([Forig fs f1], [Fpre fs f1]);
        copyfile([Forig fs f2], [Fpre fs f2]);
    end
    
end
