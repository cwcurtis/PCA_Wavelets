function nls_tester_script(Llx,K,k0,sig,tf,dt)
    
    nlvls = 6;
    dmdtol = 4;
    noise = .5;
    KT = 2*K;
        
    nsol = nls_solver_stndalne(k0,K,Llx,sig,tf,dt,noise);
    %psol = wvlt_mtlb_bltin(nsol,nlvls);
    %dmd_maker(psol,dt,dmdtol)
    
    %nmode = pca_analysis(nsol,2);
    psol = wvlt_decomp_tseries(nsol,nlvls,Llx,KT);
    dmd_maker(psol,dt,dmdtol)
    
    %dmd_maker(nsol,dt,dmdtol)