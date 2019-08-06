function nls_tester_script(Llx,K,k0,sig,tf,dt)
    
    nlvls = 3;
    KT = 2*K;
    Xmesh = linspace(-Llx,Llx,KT+1);
    Xmesh = Xmesh(1:KT)';
    osamp = 2;
    
    rphi = phirscl(vphi,Llx,KT,dx,osamp); 
    [nfin,hfltr,hflti,osamp] = nls_solver_stndalne(rphi,k0,K,Llx,sig,tf,dt);
    
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
    
    [dwtcr,dwtci,hfltrr,hfltir,gfltrr,gfltir] = wvlt_decomp_nls(nfin,gfltr,gflti,hfltr,hflti,nlvls,Llx,KT,osamp);
    
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
        ar = idwt(ar,dr,gfltrr,hfltrr);
        ai = idwt(ai,di,gfltir,hfltir);
        tot = tot*2;
        are = interp1(Xmesh(1:2^(nlvls+1-jj):KT),ar,Xmesh,'spline');
        aie = interp1(Xmesh(1:2^(nlvls+1-jj):KT),ai,Xmesh,'spline');
        figure(jj+1)
        plot(Xmesh,abs(are+1i*aie),'k-','LineWidth',2)
        
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
    