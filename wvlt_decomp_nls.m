function [dwtcr,dwtci] = wvlt_decomp_nls(f0,gfltr,gflti,hfltr,hflti,nlvls,KT)

tot = 2*KT;
dwtcr = zeros(tot,1);
dwtci = zeros(tot,1);
pcnt = KT-1;
lstp = tot-pcnt;
rstp = tot;
f0r = real(f0);
f0i = imag(f0);
dwtcr(lstp:rstp) = f0r;
dwtci(lstp:rstp) = f0i;

dwtmode('per')

for jj=1:nlvls
    [ar,dr] = dwt(f0r,gfltr,hfltr);  
    [ai,di] = dwt(f0i,gflti,hflti);  
    
    stp = KT*2^(-jj) - 1;
    rstp = lstp-1;
    lstp = rstp - stp;
    dwtcr(lstp:rstp) = dr;
    dwtci(lstp:rstp) = di;
    
    f0r = ar;
    f0i = ai;
end

dwtcr(1:lstp-1) = ar;
dwtci(1:lstp-1) = ai;
