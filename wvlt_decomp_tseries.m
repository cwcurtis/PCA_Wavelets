function psol = wvlt_decomp_tseries(nsol,nlvls,Llx,KT,osamp)

    tot = 2*KT;
    pcnt = KT-1;    
    lstp = tot-pcnt;
    
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
    
    hfunr = Kmat*hflt;
    gfunr = Kmat*gflt;

    hfltrr = real(Kmati*hfunr/KT);
    gfltrr = real(Kmati*gfunr/KT);
    
    [~,nsteps] = size(nsol);
    psol = zeros(KT,nsteps);
    
    dwtmode('per')
    
    for ll = 1:nsteps
    
        f0r = real(nsol(:,ll));
        f0i = imag(nsol(:,ll));
    
        f0r = real(ifft(fft(f0r).*conj(psif)))/sqrt(dx);
        f0i = real(ifft(fft(f0i).*conj(psif)))/sqrt(dx);
        
        for jj=1:nlvls
            [ar,~] = dwt(f0r,gfltrr,hfltrr);  
            [ai,~] = dwt(f0i,gfltrr,hfltrr);    
            stp = KT*2^(-jj) - 1;
            rstp = lstp - 1;
            lstp = rstp - stp;
            f0r = ar;
            f0i = ai;
        end

        [afr,afi] = wvlt_recon_tseries(nlvls,KT,ar,ai,gfltrr,hfltrr);

        %[rphi,~] = phirscl(afr+1i*afi,dx);
        psol(:,ll) = afr+1i*afi;
    end