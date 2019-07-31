function hflt = filter_maker(vphi,KT,dx,Llx)
    % Nkp = 2^J, where J is the number of levels in the MRA
    %{
    K = KT/2;
    odx = pi/(2*dx);
    
    rphi = phirscl(vphi,Llx,KT,dx);
    smat = repmat(Xmesh,KT,1) - 2*repmat(Xmesh',1,KT);
    hflt = zeros(KT,1);
    for nn = 1:KT
        pmat = sinc(odx*( (nn-K) + smat ) );
        hflt(nn) = sqrt(dx)*(rphi'*(pmat*rphi))/sqrt(2);
    end
    %}
    
    
    dk = pi/(KT*dx);
    K = KT/2;
    tgam = phirscl(vphi,Llx,KT,dx);
    
    Xmesh = (-Llx:dx:Llx-dx);
    Kmesh = (-pi/(2*dx):dk:pi/(2*dx)-dk)';
    Kmeshr = 
    Kmat = exp(-1i*Kmesh*Xmesh);
    dKmat = exp(-1i*2*Kmesh*Xmesh);
    hfun = sqrt(2)*((dKmat*vphi)./(Kmat*vphi)).*tgam(1:end-1);
    plot(Kmesh,abs(hfun))
    pause
    Nvals = (-K+1:K)';
    hflt = real(1/KT*exp(1i*Nvals*Kmesh')*hfun);
    
    %{
    K = KT/2;
    tgam = phirscl(vphi,Llx,KT,dx);
    
    dk = pi/(KT*dx);
    Nvals = (-K+1:K);
    Kmesh = (-pi/(2*dx):dk:pi/(2*dx))';
    %Kmeshr = (-pi:2*pi/KT:pi)';
    Xmesh = (-Llx:dx:Llx-dx);    
    Kmat1 = exp(-1i*Kmesh*Xmesh);
    Kmat2 = exp(1i*Kmesh*Xmesh/2);
    v1 = Kmat1*vphi;
    v2 = Kmat2*vphi;
    vcphi = v1.*v2.*tgam;    
    hflt = zeros(KT,1);

    for nn = 1:KT
        nscl = exp(1i*Kmesh*Nvals(nn)).*vcphi;        
        hflt(nn) = dx/sqrt(2)*real((nscl(1)+nscl(end)+2*sum(nscl(2:KT)))/(2*KT));
    end
    plot(Nvals,hflt)
    pause
    %}   
end