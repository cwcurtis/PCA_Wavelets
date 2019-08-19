function [dwtcr,dwtci,gfltrr,hfltrr,psif] = wvlt_decomp_tseries(nsol,nlvls,Llx,KT,osamp)

    tot = 2*KT;
    dwtcr = zeros(tot,1);
    dwtci = zeros(tot,1);
    pcnt = KT-1;    
    lstp = tot-pcnt;
    rstp = tot;

    K = KT/2;
    Kosmp = K*osamp;
    dx = 2*Llx/KT;
    dk = 2*pi/(dx*KT);
    Kmesh = (-pi/dx:dk:pi/dx-dk);
    Nvals = (-Kosmp+1:Kosmp);
    Nvalsr = (-K+1:K);
    Kmat = exp(-1i*Kmesh'*Nvals*dx);    
    Kmati = exp(1i*Nvalsr'*Kmesh*dx);
    
    [gflt,hflt,psif] = filter_maker(KT,dx,Llx,osamp);
    
    [~,nsteps] = size(nsol);
    
    for ll = 1:nsteps
    
    f0r = real(f0);
    f0i = imag(f0);

    
    f0r = real(ifft(fft(f0r).*conj(psif)))/sqrt(dx);
    f0i = real(ifft(fft(f0i).*conj(psif)))/sqrt(dx);
    dwtcr(lstp:rstp) = f0r;
    dwtci(lstp:rstp) = f0i;
    
    hfunr = Kmat*hflt;
    gfunr = Kmat*gflt;

    hfltrr = real(Kmati*hfunr/KT);
    gfltrr = real(Kmati*gfunr/KT);
    
    dwtmode('per')

    for jj=1:nlvls
        [ar,dr] = dwt(f0r,gfltrr,hfltrr);  
        [ai,di] = dwt(f0i,gfltrr,hfltrr);    
        stp = KT*2^(-jj) - 1;
        rstp = lstp - 1;
        lstp = rstp - stp;
        dwtcr(lstp:rstp) = dr;
        dwtci(lstp:rstp) = di;
        f0r = ar;
        f0i = ai;
    end

    dwtcr(1:lstp-1) = ar;
    dwtci(1:lstp-1) = ai;