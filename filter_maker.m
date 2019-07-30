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
    
    %{
    dk = pi/(KT*dx);
    K = KT/2;
    Xmesh = (-Llx:dx:Llx-dx);
    Kmesh = (-pi/(2*dx):dk:pi/(2*dx)-dk)';
    Kmat = exp(-1i*Kmesh*Xmesh);
    dKmat = exp(-1i*2*Kmesh*Xmesh);
    hfun = sqrt(2)*(dKmat*vphi)./(Kmat*vphi);
    plot(Kmesh,abs(hfun))
    pause
    Nvals = (-K+1:K)';
    hflt = real(1/KT*exp(1i*Nvals*Kmesh')*hfun);
    %}
    
    K = KT/2;
    tgam = phirscl(vphi,Llx,KT,dx);
    
    dk = pi/(KT*dx);
    Nvals = (-K+1:K);
    Kmesh = (-pi/(2*dx):dk:pi/(2*dx))';
    plot(Kmesh,tgam)
    pause
    Xdif = repmat(0:KT-1,KT,1) - 2*repmat((0:KT-1)',1,KT);
    Mpf = exp(-1i*pi/(2*dx)*Xdif);
    Mfn = exp(1i*pi/(2*dx)*Xdif);
    Msc = exp(1i*dk*Xdif);
    hflt = zeros(KT,1);

    for nn = 1:KT
        nscl = exp(1i*Kmesh*(Nvals(nn)+Llx)).*tgam;
        tot = Mpf*nscl(1) + Mfn*nscl(end);
        Mcc = Mpf;
        for ll=2:KT
            Mcc = nscl(ll)*Mcc.*Msc;
            tot = tot + 2*Mcc;
        end
        hflt(nn) = dx*(vphi'*(tot*vphi))/(2*sqrt(2)*KT);
    end
       
end