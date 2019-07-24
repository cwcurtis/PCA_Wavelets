function tester_script(Llx,K,k0,ep,sig,tf,dt)

    [nfin,hflt] = afm_dno_solver(K,k0,ep,Llx,sig,tf,dt);
    nlvls = -log2(ep);
    KT = 2*K;
    slvls = KT;
    Xmesh = linspace(-Llx,Llx,KT+1);
    Xmesh = Xmesh(1:KT)';
    
    dwtc = wvlt_decomp(nfin,hflt,nlvls,slvls);
    
    figure(1)
    plot(Xmesh,nfin,'k-','LineWidth',2)
    
    lstp = 1;
    rstp = KT/2^(nlvls);
    tot = KT/2^(nlvls);
        
    for jj=1:nlvls
        figure(jj+1)
        plot(1:tot,dwtc(lstp:rstp),'k-','LineWidth',2)
        tot = tot*2;
        lstp = rstp + 1;
        rstp = lstp - 1 + tot;
    end
    
    