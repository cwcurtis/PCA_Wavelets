function dmd_maker(nsol,dt)

    [~,nsteps] = size(nsol);
    
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
    fbal = bvls.*(diag(Evls).^(nsteps-1));
    naprx = sum(Mds*fbal,2);
    
    disp("Formal error is:")
    disp(norm(naprx-nsol(:,end),'inf')/norm(nsol(:,end),'inf'))
    sigtrms = log10( abs( fbal/max(abs(fbal)) ) );
    
    figure(1)
    scatter(real(evls),imag(evls))
    
    figure(2)
    scatter(real(evls), sigtrms)
    
    figure(3)
    scatter(imag(evls), sigtrms)
    
    figure(4)
    plot(abs(nsol(:,end)))
    
    figure(5)
    plot(abs(naprx))
    
    figure(6)
    surf(1:kpsz,1:mdsz,log10(abs(Mds)),'LineStyle','none')
    
    figure(7)
    surf(1:kpsz,1:kpsz,log10(abs(Mds'*Mds)),'LineStyle','none')
    