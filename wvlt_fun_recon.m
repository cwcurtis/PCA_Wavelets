function rcf = wvlt_fun_recon(ac,psif,dx,KT,lvl)
    
    tlvl = 2^lvl;
    psifr = fftshift(psif);
    psifr = ifftshift(psifr(1:tlvl:KT));
    nones = (-1).^(1:(KT/tlvl))';
    rcf = 1/(sqrt(tlvl))*real(ifft(psifr.*nones.*fft(ac)));