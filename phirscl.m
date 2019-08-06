function rphi = phirscl(vphi,Llx,KT,dx,osamp)
    
    KTT = osamp*KT;
    dk = pi/(KTT*dx);
    Kmesh = (-pi/dx:dk:pi/dx);
    Xmesh = (-Llx:dx:Llx-dx);
    Kmat1 = exp(-1i*Kmesh'*Xmesh);
    pvec = exp( 1i*angle(Kmat1*vphi) );
    Kmat2 = exp(1i*Xmesh'*Kmesh).*repmat(pvec.',KT,1);
    
    rphi = 1/(3*sqrt(dx)*KTT)*(Kmat2(:,1) + Kmat2(:,KTT+1) + 4*sum(Kmat2(:,2:2:KTT),2) + 2*sum(Kmat2(:,3:2:KTT-1),2));
    
    