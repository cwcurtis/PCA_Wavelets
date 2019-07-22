function hfun = filter_maker(vphi,Llx,Ncp,Nkp)

    dx = 2*Llx/Ncp;
    dk = 2*pi/Nkp;
    tpi = 2*pi;    
    
    M = 1/(2*dx);
    
    Xmesh = (-Llx:dx:Llx-dx)';
    Kmesh = (0:dk:tpi-dk)';
    Kmat = exp(-1i*Kmesh*Xmesh');
    Kmat2 = exp(-1i*2*Kmesh*Xmesh');
    
    gam = Kmesh;
    gam2 = 2*Kmesh;
    Svec = [(M+1/2) sin(pi*(1+1/(2*M))*(1:Ncp-1))./sin(pi*(1:Ncp-1)/M)]; 
    Skn0 = sqrt(dx)*norm(vphi);
    Sk0 = dx*sqrt(real(vphi'*toeplitz(Svec)*vphi));
    
    gam(1) = Sk0;
    gam(2:end) = Skn0;
    
    zinds = mod(gam2,tpi) == 0;
    gam2(zinds) = Sk0;
    gam2(logical(1-zinds)) = Skn0;
    
    gtot = gam./gam2;
    
    hfun = ifft(((Kmat2*vphi)./(Kmat*vphi)).*gtot);
    
end