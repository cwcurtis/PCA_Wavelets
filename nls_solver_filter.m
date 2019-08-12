function [nfin,nmode] = nls_solver_filter(k0,K,Llx,sig,tf,dt)

    nsteps = floor(tf/dt);
    KT = 2*K;
    [k0,sk,Om,ad,anl] = nls_params(k0,Llx,sig);
    Kmesh = 1i*pi/Llx*[0:K -K+1:-1]';
    Kc = floor(2*K/3);
    Kuc = KT-Kc+1;
    Kc = Kc+1;
    
    Xmesh = linspace(-Llx,Llx,KT+1);
    Xmesh = (Xmesh(1:KT))';
    ep = .5;
    sltn = sqrt(2*ad/anl)*(sech(Xmesh) + .5*sech(Xmesh-5));
    noise = fft(exp(-Xmesh.^2./(2*ep^2))/(ep*sqrt(2*pi))).*exp(1i*2*pi*rand(KT,1));
    n0 = fft(sltn)+noise*ep;
    n0(Kc:Kuc) = 0;
    
    plot(Xmesh,abs(ifft(n0)),'k-','LineWidth',2)
    pause
    
    Emh = exp(1i*dt/2*ad*Kmesh.^2);
    Em = Emh.*Emh;
    
    dts = 1i*anl*dt;
    nsol = zeros(KT,nsteps+1);
    nsol(:,1) = ifft(n0);
    
    for jj=1:nsteps
        
        a = ifft(n0);
        k1 = dts*fft((abs(a).^2).*a);
        k1(Kc:Kuc) = 0;
        
        b = ifft(Emh.*(n0  + k1/2));
        k2 = dts*fft((abs(b).^2).*b);
        k2(Kc:Kuc) = 0;
        
        c = ifft(Emh.*n0 + k2/2);
        k3 = dts*fft((abs(c).^2).*c);
        k3(Kc:Kuc) = 0;
        
        d = ifft(Em.*n0 + Emh.*k3);
        k4 = dts*fft((abs(d).^2).*d);
        k4(Kc:Kuc) = 0;
        
        n0 = Em.*(n0+k1/6) + Emh.*(k2+k3)/3 + k4/6;
        n0(Kc:Kuc) = 0;
        
        nsol(:,jj+1) = ifft(n0);
        
    end  
    nfin = ifft(n0);
    navg = mean(nsol,2);
    ndif = nfin - navg;
    [U,S,~] = svd((nsol-navg));
    SD = diag(S);
    SD = SD/max(SD);
    ikp = SD>0;
    SDr = log10(SD(ikp));
    
    ip1 = sum(ndif.*conj(U(:,1)));
    ip2 = sum(ndif.*conj(U(:,2)));
    ip3 = sum(ndif.*conj(U(:,3)));
    ip4 = sum(ndif.*conj(U(:,4)));
    
    %nmode = navg + ip1*U(:,1) + ip2*U(:,2) + ip3*U(:,3) + ip4*U(:,4);
    %nmode = navg + ip1*U(:,1);
    
    dx = Llx/K;
    nmodet = 1/(dx*sqrt(2*pi)).*exp(-Xmesh.^2/(2*dx^2));
    nmode = nmodet.*exp(1i.*Xmesh);
    figure(1)
    plot(SDr,'k-','LineWidth',2)
    figure(2)
    plot(Xmesh,abs(nmode),'k-','LineWidth',2)
    pause