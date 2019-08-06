function nls_dmd_script(Llx,K,k0,sig,tf,dt)
    
    KT = 2*K;
    osamp = 8;
    Xmesh = linspace(-Llx,Llx,KT+1);
    Xmesh = Xmesh(1:KT)';
    dx = Llx/K;
    nsteps = tf/dt;
    
    rphi = phirscl(sech(Xmesh),Llx,KT,dx,osamp); 
    
    nsol = nls_solver_stndalne(rphi,k0,K,Llx,sig,tf,dt);
    
    disp("Condition number of nsol is:")
    disp(cond(nsol(:,1:end-1)))
    
    [U,S,V] = svd(nsol(:,1:end-1),'econ');
    Sd = diag(S);
    mSv = max(Sd);
    kpinds = log10(Sd/mSv) > -18;
    Si = diag(1./Sd(kpinds));
    Amat = U(:,kpinds)'*nsol(:,2:end)*V(:,kpinds)*Si;
    
    figure(1)
    spy(U(:,kpinds)*Amat*U(:,kpinds)')
    pause
    
    figure(2)
    spy(nsol(:,1:end))
    pause
    
    [Mds,Evls] = eig(Amat);
    Mds = U(:,kpinds)*Mds;
    bvls = Mds\nsol(:,1);
    
    evls = log(diag(Evls))/dt;
    fbal = bvls.*(diag(Evls).^nsteps);
    naprx = sum(Mds*fbal,2);
    
    disp("Formal error is:")
    disp(norm(naprx-nsol(:,end))/norm(nsol(:,end)))
    
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
    plot(-K+1:K,abs(fftshift(nsol(:,end))),'r-','LineWidth',2)
    
    figure(5)
    plot(-K+1:K,abs(fftshift(naprx)),'k-','LineWidth',2)
    
    figure(6)
    plot(log10(Sd(kpinds)))
    