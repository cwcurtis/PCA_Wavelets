function nsol = nls_solver_stndalne(k0,K,Llx,sig,tf,dt)

    nsteps = floor(tf/dt);
    KT = 2*K;
    [k0,sk,Om,ad,anl] = nls_params(k0,Llx,sig);
    Kmesh = 1i*pi/Llx*[0:K -K+1:-1]';
    Kc = floor(2*K/3);
    Kuc = KT-Kc+1;
    Kc = Kc+1;
    
    Xmesh = linspace(-Llx,Llx,KT+1);
    Xmesh = (Xmesh(1:KT))';
    ep = .05;
    sltn = sqrt(2*ad/anl)*(sech(Xmesh) + ep*sech(Xmesh-5));
    noise = fft(exp(-Xmesh.^2./(2*ep^2))/(ep*sqrt(2*pi))).*exp(1i*2*pi*rand(KT,1));
    n0 = fft(sltn)+noise*ep;
    n0(Kc:Kuc) = 0;
    
    nsol = zeros(KT,nsteps+1);
    nsol(:,1) = ifft(n0);
    
    Emh = exp(1i*dt/2*ad*Kmesh.^2);
    Em = Emh.*Emh;
    
    dts = 1i*anl*dt;
    
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
    
    