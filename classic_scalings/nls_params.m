function [k0,sk,Om,ad,anl] = nls_params(k0,Llx,sig)

    k0 = k0*pi/Llx;
    
    sk = sign(k0);

    Om = sqrt(abs(k0)*(1+sig*k0^2));

    cg = (1+3*sig*k0^2)/(2*sk*Om);

    ad = (cg^2 - 3*abs(k0)*sig)/(2*Om);
    
    anl = k0/((2*sk*Om)*(4*Om^2-sk*(2*k0*(1+4*sig*k0^2))))*(sk*k0^3*(8 + sig*k0^2 + 2*(sig*k0^2)^2));
    