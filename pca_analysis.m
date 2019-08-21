function nmode = pca_analysis(nsol,tol)
    
    nfin = nsol(:,end);
    navg = mean(nsol,2);
    ndif = nfin - navg;
    [U,S,~] = svd((nsol-navg));
    SD = diag(S);
    SD = log10(SD/max(SD));
    ikp = SD>-tol;    
    ips = U(:,ikp)'*ndif;
    nmode = navg + U(:,ikp)*ips;
    
end

