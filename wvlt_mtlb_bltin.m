function psol = wvlt_mtlb_bltin(nsol,nlvls)

    [KT,nsteps] = size(nsol);
    psol = zeros(KT+2*(nlvls+1),nsteps);
    
    dwtmode('per')
    for ll=1:nsteps
        psol(1:KT,ll) = nsol(:,ll);
        [cwltsr,~] = wavedec(real(nsol(:,ll)),nlvls,'coif4');
        [cwltsi,~] = wavedec(imag(nsol(:,ll)),nlvls,'coif4');
        lstp = 1;
        rstp = KT/2^nlvls;
        cvals = cwltsr(lstp:rstp)+1i*cwltsi(lstp:rstp);
        psol(KT+1,ll) = norm(cvals,2);
        psol(KT+nlvls+1+1,ll) = norm(cvals,4);
        for mm=2:nlvls+1
            lstp = rstp + 1;
            rstp = lstp - 1 + KT/2^(nlvls-mm+2);
            cvals = cwltsr(lstp:rstp)+1i*cwltsi(lstp:rstp);
            psol(KT+mm,ll) = norm(cvals,2);
            psol(KT+nlvls+1+mm,ll) = norm(cvals,4);
        end
    end