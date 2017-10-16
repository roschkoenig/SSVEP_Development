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

Fimport     = '/Volumes/ANRI/DevEEG';



subs = cellstr(spm_select('List', Fimport, 'dir'));
textprogressbar('Copying across: ');

for s = 39:1:length(subs)
    sub = subs{s};
    textprogressbar(100*s / length(subs))
    Feeg        = [Fdata fs 'EEG' fs sub fs 'EEG' fs 'preprocessed' fs 'mat_format'];
    made        = mkdir(Feeg);
    
	tofolder    = Feeg;
    fromfolder  = [Fimport fs sub fs 'EEG' fs 'preprocessed' fs 'mat_format'];
    copyfile(fromfolder, tofolder); 
end
textprogressbar(' Done');