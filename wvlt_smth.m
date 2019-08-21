function fout = wvlt_smth(fin,nlvls)

    KT = length(fin);
    
    lstp = KT/2^(nlvls-1);
    rstp = lstp + KT*(1 - 1/(2^(nlvls-1)));
    dtlscfs = ( lstp + 1 ):rstp;
    dwtmode('per')
    [cwltsr,lv1] = wavedec(real(fin),nlvls,'coif4');
    [cwltsi,lv2] = wavedec(imag(fin),nlvls,'coif4');
    cwltsr(dtlscfs) = 0;
    cwltsi(dtlscfs) = 0;
    foutr = waverec(cwltsr,lv1,'coif4');
    fouti = waverec(cwltsi,lv2,'coif4');
    fout = foutr + 1i*fouti;
    
    