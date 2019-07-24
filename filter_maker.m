function hflt = filter_maker(vphi,Llx,Nkp)
    % Nkp = 2^J, where J is the number of levels in the MRA
    
    Ncp = length(vphi);
    dx = 2*Llx/Ncp;
    dk = 2*pi/Nkp;
    tpi = 2*pi;    
    odx = 1/(2*dx);
    M = floor(odx);
    tht = odx - M;
    
    Xmesh = (-Llx:dx:Llx-dx)';
    Kmesh = (0:dk:tpi-dk)';
    
    gvec = zeros(Nkp,1);
    for jj = 1:Nkp
        gvec(jj) = gam_comp(dx,vphi,Ncp,M,tht,Kmesh(jj));
    end
    gdvec = [gvec(1:2:end-1);gvec(1:2:end-1)];
    gtot = gdvec./gvec;
    
    Kmat = exp(-1i*Kmesh*Xmesh');
    Kmat2 = exp(-1i*2*Kmesh*Xmesh');
    
    hflt = sqrt(2)*ifft(((Kmat2*vphi)./(Kmat*vphi)).*gtot);
    
end