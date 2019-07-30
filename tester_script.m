function tester_script(Llx,K,k0,ep,sig,tf,dt)

    [nfin,hflt] = afm_dno_solver(K,k0,ep,Llx,sig,tf,dt);
    disp(sum(hflt))
    nlvls = -log2(ep);
    KT = 2*K;
    slvls = KT;
    dk = 2*pi/KT;
    dx = 2*Llx/KT;
    Xmesh = linspace(-Llx,Llx,KT+1);
    Xmesh = Xmesh(1:KT)';
    Kmesh = (-pi:dk:pi-dk);
    Nvals = (-K+1:K);
    Kmat = exp(-1i*Kmesh'*Nvals);
    
    gflt = (-1).^(-K+1:K)'.*flipud(hflt);
    gflt = [gflt(end);gflt(1:end-1)];
    hfun = Kmat*hflt;
    gfun = exp(-1i*Kmesh').*conj(exp(-1i*(Kmesh+pi)'*Nvals)*hflt);
    %gfun = fft(fftshift(gflt));
    
    gfltt = qmf(hflt,1);
    gfunt = fft((gfltt));
    %disp(norm(gflt-gfltt))
    figure(1)
    plot(1:KT,abs(hfun).^2 + abs(gfun).^2,'k-','LineWidth',2)
    pause
    %gflt = qmf(hflt,0);
    
    
    dwtc = wvlt_decomp(nfin,gflt,hflt,nlvls,slvls);
    
    figure(1)
    plot(Xmesh,nfin,'k-','LineWidth',2)
    
    lstp = 1;
    rstp = KT*ep;
    tot = KT*ep;
    figure(2)
    plot(1:tot,dwtc(lstp:rstp),'k-','LineWidth',2)
    
    lstp = rstp + 1;
    rstp = 2*rstp;
    
    for jj=2:nlvls+1
        figure(jj+1)
        plot(1:tot,dwtc(lstp:rstp),'k-','LineWidth',2)
        tot = tot*2;
        lstp = rstp + 1;
        rstp = lstp - 1 + tot;
    end
    
    a = dwtc(1:KT*ep);
    d = zeros(length(a),1);
    napprox = interpft(idwt(a,d,gflt,hflt),KT);
    
    figure(nlvls+3)
    plot(Xmesh,napprox,'k-','LineWidth',2)
    