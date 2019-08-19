function psol = wvlt_trns_proc(nsol,nlvls,Llx,KT,osamp)

    
    [dwtcr,dwtci,gfltrr,hfltrr,psif] = wvlt_decomp_tseries(nsol,nlvls,Llx,KT,osamp);
    [ar,ai] = wvlt_recon_nls(nlvls,KT,dwtcr,dwtci,gfltrr,hfltrr);
    acfr = wvlt_fun_recon(ar,psif,dx,KT,0);
    acfi = wvlt_fun_recon(ai,psif,dx,KT,0);
    af = acfr + 1i*acfi;
    [rphi,~] = phirscl(af,dx);
    