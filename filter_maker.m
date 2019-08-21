function [gflt,hflt,rphif] = filter_maker(vphi,KT,dx,Llx)
    
    dk = 2*pi/(KT*dx);
    
    %vphi = 1/(dx*sqrt(2*pi)).*exp(-Xmesh.^2/(2*dx^2));
    
    vphira = wvlt_smth(real(vphi),6);    
    vphiia = wvlt_smth(imag(vphi),6);
    avphia = abs(vphira + 1i*vphiia);
    
    [~,rphif] = phirscl(avphia,dx);
    
    Kmesh = (-pi/dx:dk:pi/dx-dk);
    mask = zeros(KT,1);
    inds = abs(Kmesh)<= pi/(2*dx);
    mask(inds) = 1;
        
    Xmesh = (-Llx:dx:Llx-dx);
    
    Kmat1 = exp(-1i*Kmesh'*Xmesh);
    Kmat2 = exp(-1i*2*Kmesh'*Xmesh);
    pvec = sqrt(2).*mask.*exp( 1i*( (Kmesh').^8.*( angle(Kmat2*avphia) - angle(Kmat1*avphia) ) ) );
    
    gflt = ifft(ifftshift(pvec));
    hflt = qmf(gflt);    
    
end