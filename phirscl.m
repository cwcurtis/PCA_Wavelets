function rphi = phirscl(vphi,KT,dx)
    dk = 2*pi/KT;
    tpi = 2*pi;    
    odx = 1/(2*dx);
    M = floor(odx);
    tht = odx - M;
    
    Kmesh = (-pi:dk:pi);
    Nvals = (0:KT-1)';
    gvec = zeros(KT+1,1);
    for jj = 1:KT+1
        kval = mod(Kmesh(jj)/dx,tpi);
        gvec(jj) = gam_comp(dx,vphi,KT,M,tht,kval);
    end
    Kmat = exp(1i*Nvals*Kmesh(2:end-1));
    rcg = 1./gvec;
    Kvec = 1/(2*KT)*( (-1).^(Nvals)*(rcg(1)+rcg(end)) + 2*Kmat*rcg(2:end-1));
    rphi = toeplitz(Kvec')*vphi;
    rphi = real(rphi); %Note, imaginary part was found to be less than machine precision.
end

