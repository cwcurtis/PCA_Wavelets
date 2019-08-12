function nls_dmd_script(Llx,K,k0,sig,tf,dt)
    
    KT = 2*K;
    Xmesh = linspace(-Llx,Llx,KT+1);
    Xmesh = Xmesh(1:KT)';
    dx = Llx/K;
    nsteps = tf/dt;
    Kc = floor(2*K/3);
    Kuc = KT-Kc+1;
    Kc = Kc+1;
    pslice = [1:Kc-1 Kuc+1:KT];
    
    %rphi = phirscl(1/(dx*sqrt(2*pi)).*exp(-Xmesh.^2/(2*dx^2)).*exp(1i*Xmesh),dx); 
    %rphi = phirscl(sech(Xmesh),dx); 
    nsol = nls_solver_stndalne(k0,K,Llx,sig,tf,dt);
    
    %disp("Condition number of nsol is:")
    %disp(cond(nsol(:,1:end-1)))
    
    [U,S,V] = svd(nsol(:,1:end-1),'econ');
    Sd = diag(S);
    mSv = max(Sd);
    kpinds = log10(Sd/mSv) > -6;
    Si = diag(1./Sd(kpinds));
    Amat = U(:,kpinds)'*nsol(:,2:end)*V(:,kpinds)*Si;
        
    [Mds,Evls] = eig(Amat);
    Mds = U(:,kpinds)*Mds;
    [mdsz,kpsz] = size(Mds);
    disp("Condition number of modes matrix is:")
    disp(cond(Mds'*Mds))
    
    bvls = Mds\nsol(:,1);
    
    evls = log(diag(Evls))/dt;
    fbal = bvls.*(diag(Evls).^nsteps);
    naprx = sum(Mds*fbal,2);
    
    disp("Formal error is:")
    %disp(norm(naprx-nsol(:,end))/norm(nsol(:,end)))
    disp(norm(naprx-nsol(:,end),'inf')/norm(nsol(:,end),'inf'))
    sigtrms = log10( abs( fbal/max(abs(fbal)) ) );
    %kptrms = sigtrms > -16;
    %disp("Number of kept terms:")
    %disp(sum(kptrms))
    
    figure(1)
    scatter(real(evls),imag(evls))
    
    figure(2)
    scatter(real(evls), sigtrms)
    
    figure(3)
    scatter(imag(evls), sigtrms)
    
    figure(4)
    plot(Kuc-KT+1:(Kc-1),abs((nsol(:,end))),'r-','LineWidth',2)
    
    figure(5)
    plot(Kuc-KT+1:(Kc-1),abs(naprx),'k-','LineWidth',2)
    
    figure(6)
    surf(1:kpsz,1:mdsz,log10(abs(Mds)),'LineStyle','none')
    
    figure(7)
    surf(1:kpsz,1:kpsz,log10(abs(Mds'*Mds)),'LineStyle','none')
    