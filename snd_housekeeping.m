function D = snd_housekeeping(e)

if strcmp('MACI64', computer)
    Fbase = '/Users/roschkoenig/Dropbox/Research/Friston Lab/1701 Sensory Network Dev';
elseif strcmp('PCWIN64', computer)
    Fbase = 'C:\Users\rrosch\Dropbox\Research\Friston Lab\1701 Sensory Network Dev';
end

if      e == 'contrast';     experiment = 'Contrast_Detection';
elseif  e == 'surround';     experiment = 'Surround_Inhibition'; 
end

fs          = filesep;
Fscripts    = [Fbase fs 'Scripts'];
Fdata       = [Fbase fs 'Data'];
Fgenex      = [Fdata fs 'Gene Expression'];
Fanalysis   = [Fbase fs 'Matlab Files' fs experiment];
Fmeeg       = [Fanalysis fs 'MEEG'];

addpath(genpath(Fscripts));

% Package up for export
%--------------------------------------------------------------------------
D.Fbase     = Fbase;
D.Fscripts  = Fscripts;
D.Fdata     = Fdata;
D.Fgenex    = Fgenex;
D.Fanalysis = Fanalysis;
D.Fmeeg     = Fmeeg;
