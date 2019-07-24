function tester_script(Llx,Ncp,k0,sig)

    Xmesh = linspace(-Llx,Llx,Ncp+1);
    Xmesh = Xmesh(1:Ncp)';
    [k0,sk,Om,ad,anl] = nls_params(k0,Llx,sig);
    n0 = sqrt(2*ad/anl)*sech(Xmesh);
    
    vphi = nls_solver(ad,anl,n0,Ncp/2,Llx,1);
    rphi = phirscl(vphi,Llx,Ncp);
    
    figure(1)
    plot(Xmesh,abs(vphi),'k-','LineWidth',2)
    figure(2)
    plot(Xmesh,abs(rphi),'k-','LineWidth',2)