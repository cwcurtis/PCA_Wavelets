% This function generates the Gamma function used throughout the algorithm.
%odx = 1/(2*dx);
%M = floor(odx);
%tht = odx - M;
    
function gam = gam_comp(dx,vphi,Ncp,M,tht,kval)

    tpi = 2*pi;
    krsc = kval/tpi;
    Nvls = dx*(1:Ncp-1);
    sNvls = dx./sin(pi*Nvls);
    
    if (krsc<tht) && (krsc+tht<1)
        Svec = exp(-1i*(kval-pi)*[0 Nvls]).*[dx*2*M sNvls.*sin(tpi*M*Nvls)];
    elseif (krsc<tht) && (krsc+tht>=1)
        Svec = exp(-1i*kval*[0 Nvls]).*[dx*2*(M+1/2) sNvls.*sin(tpi*(M+1/2)*Nvls)];
    elseif (krsc>=tht) && (krsc+tht<1)
        Svec = exp(-1i*(kval-tpi)*[0 Nvls]).*[dx*2*(M+1/2) sNvls.*sin(tpi*(M+1/2)*Nvls)];
    else
        Svec = exp(-1i*(kval-pi)*[0 Nvls]).*[dx*2*(M+1) sNvls.*sin(tpi*(M+1)*Nvls)];
    end
    
    gam = sqrt(real(vphi'*(toeplitz(Svec)*vphi)));
    
end