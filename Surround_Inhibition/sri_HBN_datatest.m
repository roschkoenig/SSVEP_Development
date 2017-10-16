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

EEGpath     = 'EEG';  
load([Fdata fs EEGpath fs 'NDARZJ460RWZ' fs 'EEG' fs 'preprocessed' fs 'mat_format' fs 'SurroundSupp_Block2.mat']);

d = result.data';
for i = 1:length(result.event);
    e(i) = str2double(result.event(i).type);
end


otrigs = [4 8];

for ot = 1:length(otrigs)

o   = [result.event(e == otrigs(ot)).sample];
 
clear allp
for k = 1:length(o)
    [p Hz] = sri_fourier(d([0:499]+o(k),:), 500);
    allp(k,:,:) = p;
end

m_trials    = squeeze(mean(allp,1));
m_chans     = mean(m_trials,2);


subplot(2,2,1 + 2*(ot-1))
plot(Hz, m_trials); hold on
title(['Onset trigger: ' num2str(otrigs(ot))]);
xlim([5 110])
xlabel('Frequency [Hz]');
ylabel('Fourier power');

subplot(2,2,[2 4])
plot(Hz, log(m_chans)); hold on
xlim([5 110])
xlabel('Frequency [Hz]');
ylabel('log power');
end

subplot(2,2,[2 4])
legend({'Trigger 4', 'Trigger 8'})

