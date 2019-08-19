function [ar,ai] = wvlt_recon_tseries(nlvls,KT,ari,aii,gfltrr,hfltrr)

    scfc = 2^nlvls;
    rstp = KT/scfc;
    tot = KT/scfc;        
    
    ar = ari;
    ai = aii;   
    
    lstp = rstp + 1;
    rstp = 2*rstp;
    
    for jj=2:nlvls+1
        ar = idwt(ar,zeros(rstp-lstp+1,1),gfltrr,hfltrr); % Filter out fast scales.
        ai = idwt(ai,zeros(rstp-lstp+1,1),gfltrr,hfltrr); % Filter out fast scales.
        
        tot = tot*2;
        
        lstp = rstp + 1;
        rstp = lstp - 1 + tot;
    end