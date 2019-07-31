function hflt = filter_maker(vphi,KT,dx,Llx)
    dk = 2*pi/(KT*dx);
    Kmesh = (-pi/dx:dk:pi/dx-dk)';
    Xmesh = (-Llx:dx:Llx-dx);
    Kmat1 = exp(-1i*Kmesh*Xmesh);
    Kmat2 = exp(-1i*Kmesh*Xmesh/dx);
    hfun = 1/sqrt(dx)*exp(1i*angle(Kmat2*vphi)).*exp(-1i*angle(Kmat1*vphi));
    kinds = abs(Kmesh) > pi;
    hfun(kinds) = 0;
    plot(Kmesh,abs(hfun))
    pause
    hflt = real(1/KT*exp(1i*dx*Kmesh*(-KT/2+1:KT/2))*hfun);
end