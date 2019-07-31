function tcg = phirscl(vphi,Llx,KT,dx)
    dk = 2*pi/KT;
    tpi = 2*pi;    
    odx = 1/(2*dx);
    M = floor(odx);
    tht = odx - M;
    
    Kmesh = (-pi:dk:pi);
    %Xmesh = linspace(-Llx,Llx,KT+1);
    %Xmesh = (Xmesh(1:KT))';
    %Nvals = (0:KT-1)';
    gvec = zeros(KT+1,1);
    dgvec = zeros(KT+1,1);
    for jj = 1:KT+1
        kval = mod(Kmesh(jj)/dx,tpi);
        dkval = mod(2*Kmesh(jj)/dx,tpi);
        gvec(jj) = gam_comp(dx,vphi,KT,M,tht,kval);
        dgvec(jj) = gam_comp(dx,vphi,KT,M,tht,dkval);
    end
    %Kmat = exp(1i*Nvals*Kmesh(2:end-1));
    tcg = gvec./dgvec;
    
    %Kvec = 1/(2*KT)*( (-1).^(Nvals)*(rcg(1)+rcg(end)) + 2*Kmat*rcg(2:end-1));
    %rphi = toeplitz(Kvec')*vphi;
    %rphi = real(rphi); %Note, imaginary part was found to be less than machine precision.
    
    %Testing
    %{
    K = KT/2;
    cplot = zeros(2*K-2,1);
    for nn=1:K-1
       rvec = sinc(pi*(nn-K)/dx + pi*(0:KT-1));
       cvec = sinc(pi*(nn-K)/dx - pi*(0:KT-1));
       rval = rphi'*(toeplitz(cvec,rvec)*rphi);
       cplot(nn) = dx*rval;
    end
    for nn=K+1:KT
       rvec = sinc(pi*(nn-K)/dx + pi*(0:KT-1));
       cvec = sinc(pi*(nn-K)/dx - pi*(0:KT-1));
       rval = rphi'*(toeplitz(cvec,rvec)*rphi);
       cplot(nn-1) = dx*rval;
    end
    disp(dx*sum(rphi))
    
    figure(1)
    plot(Xmesh,rphi,'r-','LineWidth',2)
    
    figure(2)
    plot(Kmesh,gvec,'k-','LineWidth',2)
    
    figure(3)
    plot(log10(abs(cplot)))
    
    pause
    %}
end

