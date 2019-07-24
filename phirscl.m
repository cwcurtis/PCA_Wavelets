function rphi = phirscl(vphi,Llx,Nkp)
    
    Ncp = length(vphi);
    dx = 2*Llx/Ncp;
    dk = 2*pi/Nkp;
    tpi = 2*pi;    
    odx = 1/(2*dx);
    M = floor(odx);
    tht = odx - M;
    
    Kmesh = (0:dk:tpi-dk)';
    gvec = zeros(Nkp,1);
    for jj = 1:Nkp
        kval = mod(Kmesh(jj)/dx,tpi);
        gvec(jj) = gam_comp(dx,vphi,Ncp,M,tht,kval);
    end
    
    Kvec = ifft(1./gvec);
    rphi = toeplitz(Kvec')*vphi;
    disp(max(abs(imag(rphi))))
    rphi = real(rphi);
end

