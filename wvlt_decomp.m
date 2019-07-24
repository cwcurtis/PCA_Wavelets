function dwtc = wvlt_decomp(f0,hflt,nlvls,slvls)

gflt = qmf(hflt,0);

tot = slvls*(2-2^(-nlvls));
dwtc = zeros(tot,1);
pcnt = slvls-1;
lstp = tot-pcnt;
rstp = tot;
dwtc(lstp:rstp) = f0;
dwtmode('per')

for jj=1:nlvls
    [a,d] = dwt(f0,gflt,hflt);  
    %hflt = hflt(1:2:end-1);
    %gflt = qmf(hflt);
    
    stp = slvls*2^(-jj) - 1;
    rstp = lstp-1;
    lstp = rstp - stp;
    dwtc(lstp:rstp) = d;
    f0 = a;
end