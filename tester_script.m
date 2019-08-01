function tester_script(Llx,K,k0,ep,sig,tf,dt)

    [nfin,hflt] = afm_dno_solver(K,k0,ep,Llx,sig,tf,dt);
    
    disp(sum(hflt))
    nlvls = -log2(ep);
    KT = 2*K;
    slvls = KT;
    gflt = (-1).^(-K+1:K)'.*flipud(hflt);
    gflt = [gflt(KT);gflt(1:KT-1)];
    
    Xmesh = linspace(-Llx,Llx,KT+1);
    Xmesh = Xmesh(1:KT)';
    
    dx = 2*Llx/KT;
    dk = 2*pi/(dx*KT);
    Kmesh = (-pi/dx:dk:pi/dx-dk);
    Nvals = (-K+1:K);
    Kmat = exp(-1i*Kmesh'*Nvals*dx);
    hfun = Kmat*hflt;
    %gfun = exp(-1i*dx*Kmesh').*conj(exp(-1i*(Kmesh+pi/dx)'*Nvals*dx)*hflt);
    gfun = Kmat*gflt;
    plot(Kmesh,abs(hfun).^2+abs(gfun).^2)
    pause
    
    dwtc = wvlt_decomp(nfin,gflt,hflt,nlvls,slvls);
    
    figure(1)
    plot(Xmesh,nfin,'k-','LineWidth',2)
    
    lstp = 1;
    rstp = KT*ep;
    tot = KT*ep;
    
    figure(2)
    plot(Xmesh,interpft(dwtc(lstp:rstp),KT),'k-','LineWidth',2)
    
    a = dwtc(1:KT*ep);
       
    lstp = rstp + 1;
    rstp = 2*rstp;
    
    for jj=2:nlvls+1
        d = dwtc(lstp:rstp);
        a = idwt(a,d,gflt,hflt);
        tot = tot*2;
        
        figure(jj+1)
        plot(Xmesh,interpft(a,KT),'k-','LineWidth',2)
        
        lstp = rstp + 1;
        rstp = lstp - 1 + tot;
    end
    
    figure(nlvls+3)
    plot(Xmesh,nfin,'k-',Xmesh,a,'r-','LineWidth',2)
        
    figure(nlvls+4)
    plot(Xmesh,log10(abs(nfin-a)),'k-','LineWidth',2)
    
    disp("Final Error in Reconstruction is:")
    disp(norm(nfin - a)/norm(nfin))
    