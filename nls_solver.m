function dmode = nls_solver(ad,anl,n0,K,Llx,tf)

    dt = .01;
    
    nsteps = floor(tf/dt);
    nsol = zeros(2*K,nsteps+1);
    
    Kmesh = 1i*pi/Llx*[0:K -K+1:-1]';
    Linv = (1 - 1i*3*dt/4*ad*Kmesh.^2).^(-1);

    nnm1 = n0;
    nsol(:,1) = n0;
    nphys = ifft(n0);
    nln = 1i*anl*fft(nphys.^2.*conj(nphys));
    
    nlnm1 = nln;
    nlnm2 = nlnm1;
    nlnm3 = nlnm2;
    
    for jj=1:nsteps
        nphys = ifft(n0);
        nln = 1i*anl*fft(nphys.^2.*conj(nphys));
        nlstp = 55/24*nln - 59/24*nlnm1 + 37/24*nlnm2 - 3/8*nlnm3;    
        np1 = Linv.*(n0 + nnm1/3 + dt*nlstp) - nnm1/3;
        nnm1 = n0;
        n0 = np1;
        nsol(:,jj+1) = ifft(n0);
    end
    
    [U,~,~] = svd(nsol);
    dmode = U(:,1);