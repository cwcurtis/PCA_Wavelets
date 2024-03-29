function psol = wvlt_decomp_tseries(nsol,nlvls,Llx,KT)

    K = KT/2;
    dx = 2*Llx/KT;
    dk = 2*pi/(dx*KT);
    Kmesh = (-pi/dx:dk:pi/dx-dk);
    Nvals = (-K+1:K);
    Kmat = exp(-1i*Kmesh'*Nvals*dx);    
    Kmati = exp(1i*Nvals'*Kmesh*dx);
    
    [gflt,hflt,psif] = filter_maker(nsol(:,end),KT,dx,Llx);
    
    hfunr = Kmat*hflt;
    gfunr = Kmat*gflt;

    hfltrr = real(Kmati*hfunr/KT);
    gfltrr = real(Kmati*gfunr/KT);
    
    [~,nsteps] = size(nsol);
    psol = zeros(2*KT,nsteps);
    dwtcr = zeros(KT,1);
    dwtci = zeros(KT,1);
    
    dwtmode('per')
    
    for ll = 1:nsteps
    
        f0r = real(nsol(:,ll));
        f0i = imag(nsol(:,ll));
    
        f0r = real(ifft(fft(f0r).*conj(psif)))/sqrt(dx);
        f0i = real(ifft(fft(f0i).*conj(psif)))/sqrt(dx);
    
        lstp = KT+1;
               
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
        psol(:,ll) = [nsol(:,ll);dwtcr+1i*dwtci];
    end