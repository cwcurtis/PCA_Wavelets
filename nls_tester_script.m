function nls_tester_script(Llx,K,k0,sig,tf,dt)
    
    nlvls = 3;
    KT = 2*K;
    Xmesh = linspace(-Llx,Llx,KT+1);
    Xmesh = Xmesh(1:KT)';
    osamp = 4;
    dx = Llx/K;
    
    [nfin,nmode] = nls_solver_filter(k0,K,Llx,sig,tf,dt);
    hfltr = filter_maker(real(nmode),KT,dx,Llx,osamp);
    hflti = filter_maker(imag(nmode),KT,dx,Llx,osamp);
    
    KTT = KT*osamp;
    Kosmp = K*osamp;
    nonesosmp = (-1).^(-Kosmp+1:Kosmp);
    gfltr = nonesosmp'.*flipud(hfltr);
    gfltr = [gfltr(KTT);gfltr(1:KT*osamp-1)];
    gflti = nonesosmp'.*flipud(hflti);
    gflti = [gflti(KTT);gflti(1:KTT-1)];
        
    %dx = 2*Llx/KT;
    %dk = 2*pi/(dx*KT);
    %Kmesh = (-pi/dx:dk:pi/dx-dk);
    %Nvalsosmp = (-Kosmp+1:Kosmp);
    %Kmat = exp(-1i*Kmesh'*Nvalsosmp*dx);
    %hfunr = Kmat*hfltr;
    %gfunr = Kmat*gfltr;
    %hfuni = Kmat*hflti;
    %gfuni = Kmat*gflti;
    
    [dwtcr,dwtci,hfltrr,hfltir,gfltrr,gfltir,psirf,psiif] = wvlt_decomp_nls(nfin,nmode,gfltr,gflti,hfltr,hflti,nlvls,Llx,KT,osamp);
    
    figure(1)
    plot(Xmesh,abs(nfin),'k-','LineWidth',2)
    
    lstp = 1;
    scfc = 2^nlvls;
    rstp = KT/scfc;
    tot = KT/scfc;        
    
    frncr = wvlt_fun_recon(dwtcr(lstp:rstp),psirf,dx,KT,nlvls);
    frnci = wvlt_fun_recon(dwtci(lstp:rstp),psiif,dx,KT,nlvls);
    
    figure(2)
    plot(Xmesh,abs(interpft(frncr+1i*frnci,KT)),'k-','LineWidth',2)
    
    ar = dwtcr(1:KT/scfc);
    ai = dwtci(1:KT/scfc);   
    
    lstp = rstp + 1;
    rstp = 2*rstp;
    
    for jj=2:nlvls+1
        dr = dwtcr(lstp:rstp);
        di = dwtci(lstp:rstp);
        ar = idwt(ar,dr,gfltrr,hfltrr);
        ai = idwt(ai,di,gfltir,hfltir);
        frncr = wvlt_fun_recon(ar,psirf,dx,KT,nlvls-(jj-1));
        frnci = wvlt_fun_recon(ai,psiif,dx,KT,nlvls-(jj-1));
    
        tot = tot*2;
        are = interp1(Xmesh(1:2^(nlvls+1-jj):KT),frncr,Xmesh,'spline');
        aie = interp1(Xmesh(1:2^(nlvls+1-jj):KT),frnci,Xmesh,'spline');
        figure(jj+1)
        plot(Xmesh,abs(are+1i*aie),'k-','LineWidth',2)
        
        lstp = rstp + 1;
        rstp = lstp - 1 + tot;
    end
    
    af = frncr + 1i*frnci;
    
    figure(nlvls+3)
    plot(Xmesh,real(nfin),'k-',Xmesh,frncr,'r-','LineWidth',2)
        
    figure(nlvls+4)
    plot(Xmesh,imag(nfin),'k-',Xmesh,frnci,'r-','LineWidth',2)
    
    figure(nlvls+5)
    plot(Xmesh,log10(abs(nfin-af)),'k-','LineWidth',2)
    
    disp("Final Error in Reconstruction is:")
    disp(norm(nfin - af)/norm(nfin))
    