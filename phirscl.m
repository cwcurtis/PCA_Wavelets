function rphi = phirscl(dx,vphi,Ncp)
    
    gam = gam_comp(dx,vphi,Ncp);
    Kvec = ifft(1./gam);
    rphi = toeplitz(Kvec')*vphi;
    
end

