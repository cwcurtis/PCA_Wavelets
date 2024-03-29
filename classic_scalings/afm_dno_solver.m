function [nfin,hflt] = afm_dno_solver(K,k0,ep,Llx,sig,tf,dt)
    % K - number of modes in spatial Fourier series.  Typically K=128 will
    % do.
    
    % k0 - central wave number of packet.  Typically 10 will do.  
    
    % ep - magnitude of perturbations away from linearity. Typically .1
    % will do.  
    
    % Llx - domain size.  Typically Llx = 100 will do.  
    
    % sig - surface tension = 1e-5
    
    % tf - length of time to run simulation
    
    % dt - time step for simulation.  Typically .05 will do.
    
    KT = 2*K;
    % Find the wave numbers to implement the 2/3 de-aliasing throughout
    Kc = floor(2*K/3);
    Kuc = KT-Kc+1;
    Kc = Kc+1;
    
    Xmesh = linspace(-Llx,Llx,KT+1);
    Xmesh = Xmesh(1:KT)';    
    
    Kmesh = pi/Llx*[0:K -K+1:-1]';
        
    nmax = round(tf/dt);
    [k0,sk,Om,ad,anl] = nls_params(k0,Llx,sig);

    disp('Dispersion coefficient is')
    disp(ad)
    
    disp('Nonlinearity coefficient is')
    disp(anl)
    
    disp('Soliton amplitude is')
    disp(sqrt(2*ad/anl))
    
    %L1 = Kmesh.*tanh(mu.*Kmesh)/mu;
    L1 = abs(Kmesh).*(1+sig.*Kmesh.^2);
    
    Linvd = (1 + 9*dt^2/16*L1).^(-1);
    Linv12 = 3*dt/4*abs(Kmesh).*Linvd;
    Linv21 = -3*dt/4*(1+sig.*Kmesh.^2).*Linvd;
    
    % Soliton ICs
    c = 2*sqrt(ad)*sqrt(2*ad-ad);
    etan = sqrt(2*ad/anl)*sech(ep*Xmesh).*exp(1i*(k0*Xmesh - c*ep*Xmesh/(2*ad) ));
    qn = -1i*sk*Om/k0*etan;
    stdv = .01;
    
    noise = sqrt(Llx/(3*stdv*pi))*(1/(2*pi))^(.25)*( exp(-(Kmesh-k0/2).^2/(4.*stdv)).*exp(2*pi*1i*rand(KT,1)) ...
            + exp(-(Kmesh-k0).^2/(4.*stdv)).*exp(2*pi*1i*rand(KT,1)) ...
            + exp(-(Kmesh-2*k0).^2/(4.*stdv)).*exp(2*pi*1i*rand(KT,1)) );
          
    noise = 2*sqrt(2*ad/anl)*real(ifft(noise));
    etanp = 2*real(etan);
    qnp = 2*real(qn);
    etan = fft(etanp);
    qn = fft(qnp);
    
    n0 = .5*(etanp + 1i*abs(k0)*qnp/Om).*exp(-1i*k0*Xmesh);
    dmode = nls_solver(ad,anl,fft(interpft(n0,2*K*ep)),K*ep,Llx*ep,tf*ep^2);
    dxf = (2*Llx/KT);
    %fdmode = fst_samp(dmode,KT,ep);
    %fdmode = 2*real(sech(Xmesh).*fdmode.*exp(1i*k0*Xmesh/ep));
    %fdmode = .5*(fdmode + flipud(fdmode));
    fdmode = 2*real(sech(Xmesh).*exp(1i*k0*Xmesh));
    hflt = filter_maker(fdmode,KT,dxf,Llx);
    
    % Generic ICs
    % etan = fft(cos(k0*Xmesh).*sech(ep*Xmesh));
    % qn = fft(sin(k0*Xmesh).*sech(ep*Xmesh));
    
    etan(Kc:Kuc) = 0;
    qn(Kc:Kuc) = 0;
    G0 = real(ifft(L1.*qn));
    
    etanm1 = etan;
    qnm1 = qn;
    
    nln = nonlinearity(K,etan,qn,G0,ep,sig,Kmesh);
    
    nlnm1 = nln;
    nlnm2 = nlnm1;
    nlnm3 = nlnm2;
    
    for jj=1:nmax
        
        G0 = real(ifft(L1.*qn));
        nln = nonlinearity(K,etan,qn,G0,ep,sig,Kmesh);
        
        nlvecn = 55/24*nln(1:KT) - 59/24*nlnm1(1:KT) + 37/24*nlnm2(1:KT) - 3/8*nlnm3(1:KT);    
        nlvecq = 55/24*nln(KT+1:2*KT) - 59/24*nlnm1(KT+1:2*KT) + 37/24*nlnm2(KT+1:2*KT) - 3/8*nlnm3(KT+1:2*KT);
        
        eta1 = Linvd.*(etan + 1/3*etanm1 + dt*nlvecn);
        eta2 = Linv12.*(qn + 1/3*qnm1 + dt*nlvecq);
        
        q1 = Linvd.*(qn + 1/3*qnm1 + dt*nlvecq);
        q2 = Linv21.*(etan + 1/3*etanm1 + dt*nlvecn);
        
        etanp1 = -etanm1/3 + eta1 + eta2;
        qnp1 = -qnm1/3 + q1 + q2;        
        
        etanm1 = etan;
        etan = etanp1;
        
        qnm1 = qn;
        qn = qnp1;
        
        nlnm3 = nlnm2;
        nlnm2 = nlnm1;
        nlnm1 = nln;        
    end
    
    nfin = real(ifft(etan));