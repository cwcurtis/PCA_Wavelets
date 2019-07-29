function hflt = filter_maker(vphi,KT,dx,Llx)
    % Nkp = 2^J, where J is the number of levels in the MRA
    K = KT/2;
    odx = pi/(2*dx);
    
    Xmesh = (-Llx:dx:Llx-dx);
    rphi = phirscl(vphi,KT,dx);
    smat = repmat(Xmesh,KT,1) - 2*repmat(Xmesh',1,KT);
    hflt = zeros(KT,1);
    for nn = 1:KT
        pmat = sinc(odx*( (nn-K) + smat ) );
        hflt(nn) = sqrt(dx)*(rphi'*(pmat*rphi))/sqrt(2);
    end
    
end