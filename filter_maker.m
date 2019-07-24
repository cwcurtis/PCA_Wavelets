function hflt = filter_maker(vphi,Llx,K)
    % Nkp = 2^J, where J is the number of levels in the MRA
    KT = 2*K;
    Ncp = length(vphi);
    dx = 2*Llx/Ncp;
    odx = pi/(2*dx);
    
    Xmesh = (-Llx:dx:Llx-dx);
    rphi = phirscl(vphi,Llx,Ncp);
    smat = repmat(Xmesh,Ncp,1) - 2*repmat(Xmesh',1,Ncp);
    hflt = zeros(KT,1);
    
    for nn = 1:KT
        pmat = sinc(odx*( (nn-K) + smat ) );
        hflt(nn) = rphi'*(pmat*rphi)/sqrt(2);
    end
    
end