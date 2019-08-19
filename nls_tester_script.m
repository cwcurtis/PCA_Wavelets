function nls_tester_script(Llx,K,k0,sig,tf,dt)
    
    nlvls = 5;
    KT = 2*K;
    osamp = 1;
    
    nsol = nls_solver_stndalne(k0,K,Llx,sig,tf,dt);
    psol = wvlt_decomp_tseries(nsol,nlvls,Llx,KT,osamp);

    dmd_maker(psol,dt)
    