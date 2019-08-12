function fwc = wvlt_init_cond_trans(f0,psif)

    fwc = real(ifft(fft(f0).*conj(psif)));

    