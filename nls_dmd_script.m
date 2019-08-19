function nls_dmd_script(Llx,K,k0,sig,tf,dt)
    
    nsol = nls_solver_stndalne(k0,K,Llx,sig,tf,dt);
    dmd_maker(nsol,dt)