function [nfin,hfltr,hflti] = nls_solver_stndalne(k0,K,Llx,sig,tf,dt)

    KT = 2*K;
    dx = Llx/K;
    nsteps = floor(tf/dt);
    nsol = zeros(KT,nsteps+1);
    
    [k0,sk,Om,ad,anl] = nls_params(k0,Llx,sig);
    Kmesh = pi/Llx*[0:K -K+1:-1]';
    Xmesh = linspace(-Llx,Llx,KT+1);
    Xmesh = Xmesh(1:KT)';
    
    c = 2*sqrt(ad)*sqrt(2*ad-ad);
    sltn = sqrt(2*ad/anl)*sech(Xmesh).*exp(-1i*c*Xmesh/(2*ad));
    stdv = .01;
    
    noise = sqrt(Llx/(3*stdv*pi))*(1/(2*pi))^(.25)*( exp(-(Kmesh-1/2).^2/(4.*stdv)).*exp(2*pi*1i*rand(KT,1)) ...
            + exp(-(Kmesh).^2/(4.*stdv)).*exp(2*pi*1i*rand(KT,1)) ...
            + exp(-(Kmesh-2).^2/(4.*stdv)).*exp(2*pi*1i*rand(KT,1)) );
          
    noise = sqrt(2*ad/anl)*noise;
    n0 = fft(sltn);
    
    Linv = (1 + 1i*3*dt/4*ad*Kmesh.^2).^(-1);

    nnm1 = n0;
    nsol(:,1) = ifft(n0);
    nphys = ifft(n0);
    nln = 1i*anl*fft(nphys.^2.*conj(nphys));
    
    nlnm1 = nln;
    nlnm2 = nlnm1;
    nlnm3 = nlnm2;
    
    for jj=1:nsteps
        nphys = ifft(n0);
        nln = 1i*anl*fft((nphys.^2).*conj(nphys));
        nlstp = 55/24*nln - 59/24*nlnm1 + 37/24*nlnm2 - 3/8*nlnm3;    
        np1 = Linv.*(n0 + nnm1/3 + dt*nlstp) - nnm1/3;
        nnm1 = n0;
        n0 = np1;
        nsol(:,jj+1) = ifft(n0);
    end
    
    [U,S,~] = svd(nsol);
    Sd = diag(S);
    mSd = max(Sd);
    nfin = nsol(:,end);
    dmode = U(:,1);
    dmoder = real(dmode);
    %dmoder = .5*(dmoder + flipud(dmoder));
    dmodei = imag(dmode);
    %dmodei = .5*(dmodei + flipud(dmodei));
    
    hfltr = filter_maker(dmoder,KT,dx,Llx);
    hflti = filter_maker(dmodei,KT,dx,Llx);
    