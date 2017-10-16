function [P1 Hz] = sri_fourier(dat, Fs)
    fdat = fft(dat);
    L    = size(dat,1) - rem(size(dat,1), 2);
    P2   = abs(fdat / L);               % Two sided power spectrum
    P1   = P2(1:L/2+1,:);               % One sided power spectrum
    
    P1(2:end-1,:) = 2*P1(2:end-1,:);
    Hz = Fs*(0:(L/2))/L; 
end