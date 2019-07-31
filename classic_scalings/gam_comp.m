% This function generates the Gamma function used throughout the algorithm.
%odx = 1/(2*dx);
%M = floor(odx);
%tht = odx - M;
    
function gam = gam_comp(dx,vphi,KT,M,tht,kval)

    tpi = 2*pi;
    krsc = kval/tpi;
    sk = sign(krsc);
    krtt = abs(krsc) + tht;
    Nvls = dx*(1:KT-1);
    sNvls = dx./sin(pi*Nvls);
    if (tht<abs(krsc)) && (krtt<1)
        Svec = exp(-1i*(kval-sk*pi)*[0 Nvls]).*[dx*2*M sNvls.*sin(tpi*M*Nvls)];
    elseif (tht<abs(krsc)) && (krtt>=1)
        Svec = exp(-1i*kval*[0 Nvls]).*[dx*2*(M+1/2) sNvls.*sin(tpi*(M+1/2)*Nvls)];
    elseif (tht>=abs(krsc)) && (krtt<1)
        Svec = exp(-1i*(kval-sk*tpi)*[0 Nvls]).*[dx*2*(M+1/2) sNvls.*sin(tpi*(M+1/2)*Nvls)];
    elseif (tht>=abs(krsc)) && (krtt>=1)
        Svec = exp(-1i*(kval-sk*pi)*[0 Nvls]).*[dx*2*(M+1) sNvls.*sin(tpi*(M+1)*Nvls)];
    end
    tval = vphi'*(toeplitz(Svec)*vphi);
    gam = sqrt(dx*real(tval));
    
end