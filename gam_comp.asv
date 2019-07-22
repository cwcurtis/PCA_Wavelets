% This function generates the Gamma function used throughout the algorithm.

function gam = gam_comp(dx,vphi,Ncp)

    tpi = 2*pi;    
    M = 1/(2*dx);
    gam = linspace(0,tpi,Ncp+1);
    gam = gam(1:end-1)';
    gam2 = 2*gam;
    
    Svec = [(M+1/2) sin(pi*(1+1/(2*M))*(1:Ncp-1))./sin(pi*(1:Ncp-1)/M)]; 
    Skn0 = dx*sqrt(2*M)*norm(vphi);
    Sk0 = dx*sqrt(real(vphi'*toeplitz(Svec)*vphi));
    
    gam(1) = Sk0;
    gam(2:end) = Skn0;
    
    gam2(1) = Sk0;
    gam2()
    gam2(logical(1-zinds)) = Skn0;
    
end