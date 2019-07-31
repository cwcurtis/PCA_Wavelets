function dwtc = wvlt_decomp(f0,gflt,hflt,nlvls,KT)

tot = 2*KT;
dwtc = zeros(tot,1);
pcnt = KT-1;
lstp = tot-pcnt;
rstp = tot;
dwtc(lstp:rstp) = f0;
dwtmode('per')

for jj=1:nlvls
    [a,d] = dwt(f0,gflt,hflt);  
    stp = KT*2^(-jj) - 1;
    rstp = lstp-1;
    lstp = rstp - stp;
    dwtc(lstp:rstp) = d;
    f0 = a;
end

dwtc(1:lstp-1) = a;