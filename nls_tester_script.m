function nls_tester_script(Llx,K,k0,sig,tf,dt)

    [nfin,hfltr,hflti] = nls_solver_stndalne(k0,K,Llx,sig,tf,dt);
    
    nlvls = 2;
    KT = 2*K;
    slvls = KT;
    
    gfltr = (-1).^(-K+1:K)'.*flipud(hfltr);
    gfltr = [gfltr(KT);gfltr(1:KT-1)];
    gflti = (-1).^(-K+1:K)'.*flipud(hflti);
    gflti = [gflti(KT);gflti(1:KT-1)];
        
    Xmesh = linspace(-Llx,Llx,KT+1);
    Xmesh = Xmesh(1:KT)';
    
    dx = 2*Llx/KT;
    dk = 2*pi/(dx*KT);
    Kmesh = (-pi/dx:dk:pi/dx-dk);
    Nvals = (-K+1:K);
    Kmat = exp(-1i*Kmesh'*Nvals*dx);
    hfunr = Kmat*hfltr;
    gfunr = Kmat*gfltr;
    hfuni = Kmat*hflti;
    gfuni = Kmat*gflti;
    
    plot(Kmesh,abs(hfunr).^2+abs(gfunr).^2)
    pause
    
    plot(Kmesh,abs(hfuni).^2+abs(gfuni).^2)
    pause
    
    [dwtcr,dwtci] = wvlt_decomp_nls(nfin,gfltr,gflti,hfltr,hflti,nlvls,slvls);
    
    figure(1)
    plot(Xmesh,abs(nfin),'k-','LineWidth',2)
    
    lstp = 1;
    scfc = 2^nlvls;
    rstp = KT/scfc;
    tot = KT/scfc;
    
    figure(2)
    plot(Xmesh,abs(interpft(dwtcr(lstp:rstp)+1i*dwtci(lstp:rstp),KT)),'k-','LineWidth',2)
    
    ar = dwtcr(1:KT/scfc);
    ai = dwtci(1:KT/scfc);   
    
    lstp = rstp + 1;
    rstp = 2*rstp;
    
    for jj=2:nlvls+1
        dr = dwtcr(lstp:rstp);
        di = dwtci(lstp:rstp);
        ar = idwt(ar,dr,gfltr,hfltr);
        ai = idwt(ai,di,gflti,hflti);
        tot = tot*2;
        
        figure(jj+1)
        plot(Xmesh,abs(interpft(ar+1i*ai,KT)),'k-','LineWidth',2)
        
        lstp = rstp + 1;
        rstp = lstp - 1 + tot;
    end
    
    af = ar + 1i*ai;
    
    figure(nlvls+3)
    plot(Xmesh,real(nfin),'k-',Xmesh,real(af),'r-','LineWidth',2)
        
    figure(nlvls+4)
    plot(Xmesh,imag(nfin),'k-',Xmesh,imag(af),'r-','LineWidth',2)
    
    figure(nlvls+5)
    plot(Xmesh,log10(abs(nfin-af)),'k-','LineWidth',2)
    
    disp("Final Error in Reconstruction is:")
    disp(norm(nfin - af)/norm(nfin))
    