function nls_tester_script(Llx,K,k0,sig,tf,dt)
    
    nlvls = 7;
    dmdtol = 3;
    noise = .5;
    KT = 2*K;
    dx = Llx/K;
    
    nsol = nls_solver_stndalne(k0,K,Llx,sig,tf,dt,noise);
    
    [evalsnrm,strmsnrm] = dmd_maker(nsol,dt,dmdtol);
    
    psolmtlb = wvlt_mtlb_bltin(nsol,nlvls);
    [evalsmtlb,strmsmtlb] = dmd_maker(psolmtlb,dt,dmdtol);
    
    %psolemmd = wvlt_decomp_tseries(nsol,nlvls,Llx,KT);
    %[evalsemmd,strmsemmd] = dmd_maker(psolemmd,dt,dmdtol);
    
    figure(1)
    hold on 
    scatter(real(evalsnrm),imag(evalsnrm),20,'r')
    scatter(real(evalsmtlb),imag(evalsmtlb),10,'b','filled')
    %scatter(real(evalsemmd),imag(evalsemmd),10,'g','filled')
    hold off
    
    figure(2)
    hold on 
    scatter(real(evalsnrm),strmsnrm,20,'r')
    scatter(real(evalsmtlb),strmsmtlb,10,'b','filled')
    %scatter(real(evalsemmd),strmsemmd,10,'g','filled')
    hold off
    
    figure(3)
    hold on 
    scatter(imag(evalsnrm),strmsnrm,20,'r')
    scatter(imag(evalsmtlb),strmsmtlb,10,'b','filled')
    %scatter(imag(evalsemmd),strmsemmd,10,'g','filled')
    hold off
    