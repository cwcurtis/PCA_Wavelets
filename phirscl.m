function rphi = phirscl(vphi,KT,dx)
    dk = 2*pi/KT;
    tpi = 2*pi;    
    odx = 1/(2*dx);
    M = floor(odx);
    tht = odx - M;
    
    Kmesh = (-pi:dk:pi);
    %Xmesh = linspace(-Llx,Llx,KT+1);
    %Xmesh = (Xmesh(1:KT))';
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
    
    %Testing
    K = KT/2;
    nones = (-1).^(0:KT-1);
    for nn=1:K-1
       rvec = nones./(nn-K+dx*(0:KT-1));
       cvec = nones./(nn-K-dx*(0:KT-1));
       check = dx^2/pi*sin(pi*(nn-K)/dx)*rphi'*(toeplitz(cvec,rvec)*rphi);
       disp(check)
       if(isnan(check))
           disp(nn)
       end
    end
    
    %for nn=K+1:KT
    %   rvec = nones./(nn-K+dx*(0:KT-1));
    %   cvec = nones./(nn-K-dx*(0:KT-1));
    %   check = dx^2/pi*sin(pi*(nn-K)/dx)*rphi'*(toeplitz(cvec,rvec)*rphi);
    %   disp(check)
    %end
    pause
end

