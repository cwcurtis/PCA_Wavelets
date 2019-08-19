function rcf = wvlt_fun_recon(ac,psif,dx,KT,lvl)
    
    tlvl = 2^lvl;
    tpi = 2*pi;
    dk = tpi/KT;
    Ktt = KT/tlvl;
    kmesh = (0:dk:tpi-dk);
    acf = exp(-1i*kmesh'*(-Ktt/2+1:Ktt/2))*ac;
    nones = (-1).^([0:KT/2 -KT/2+1:-1]');
    rc = 1/(sqrt(dx)*sqrt(tlvl))*real(ifft(nones.*psif.*acf));
    rcf = rc(1:tlvl:KT);