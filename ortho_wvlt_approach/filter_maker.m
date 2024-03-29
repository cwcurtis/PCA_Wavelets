function [gflt,hflt,rphif] = filter_maker(KT,dx,Llx,osamp)
    
    KT = osamp*KT;
    dk = pi/(KT*dx);
    Xmesh = (-Llx:dx:Llx-dx)';
    vphi = 1/(dx*sqrt(2*pi)).*exp(-Xmesh.^2/(2*dx^2));
    [~,rphif] = phirscl(vphi,dx);
    
    Kmesh = (-pi/(2*dx):dk:pi/(2*dx));
    Xmesh = (-Llx:dx:Llx-dx);
    Kmat1 = exp(-1i*Kmesh'*Xmesh);
    Kmat2 = exp(-1i*2*Kmesh'*Xmesh);
    pvec = exp( 1i*( angle(Kmat2*vphi) - angle(Kmat1*vphi) ) );
    Kvals = 1i*dx*(-KT/2+1:KT/2)';
    
    odd = 4*exp(Kvals*Kmesh(2:2:KT))*pvec(2:2:KT);
    evn = 2*exp(Kvals*Kmesh(3:2:KT-1))*pvec(3:2:KT-1);
    gflt = real(1/(3*sqrt(2)*KT)*( exp(Kvals*Kmesh(1))*pvec(1) +  exp(Kvals*Kmesh(KT+1))*pvec(KT+1) + odd + evn ));  
    
    hflt = qmf(gflt);
    
end