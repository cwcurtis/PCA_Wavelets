function psol = wvlt_mtlb_bltin(nsol,nlvls)

    [KT,nsteps] = size(nsol);
    psol = zeros(2*KT,nsteps);
    
    lstp = KT/2^(nlvls-1);
    rstp = lstp + KT*(1 - 1/(2^(nlvls-1)));
    dtlscfs = ( lstp + 1 ):rstp;
    dwtmode('per')
    for ll=1:nsteps
       
        [cwltsr,~] = wavedec(real(nsol(:,ll)),nlvls,'coif4');
        [cwltsi,~] = wavedec(imag(nsol(:,ll)),nlvls,'coif4');
        %cwlts(dtlscfs) = 0;
        psol(:,ll) = [nsol(:,ll);cwltsr+1i*cwltsi];
        
    end