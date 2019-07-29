function fdmode = fst_samp(dmode,KT,ep)
dfreq = fftshift(fft(dmode));
fdmode = zeros(KT,1);
fdmode(1:1/ep:KT) = dfreq;
fdmode = sqrt(ep)*ifft(ifftshift(fdmode));

    