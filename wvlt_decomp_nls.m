function [dwtcr,dwtci,hfltrr,hfltir,gfltrr,gfltir,psirf,psiif] = wvlt_decomp_nls(f0,nmode,gfltr,gflti,hfltr,hflti,nlvls,Llx,KT,osamp)

tot = 2*KT;
dwtcr = zeros(tot,1);
dwtci = zeros(tot,1);
pcnt = KT-1;
lstp = tot-pcnt;
rstp = tot;

f0r = real(f0);
f0i = imag(f0);

K = KT/2;
Kosmp = K*osamp;
dx = 2*Llx/KT;
dk = 2*pi/(dx*KT);
Kmesh = (-pi/dx:dk:pi/dx-dk);
Nvals = (-Kosmp+1:Kosmp);
Nvalsr = (-K+1:K);
Kmat = exp(-1i*Kmesh'*Nvals*dx);    
Kmati = exp(1i*Nvalsr'*Kmesh*dx);

[~,psirf] = phirscl(real(nmode),dx);
[~,psiif] = phirscl(imag(nmode),dx);

fwcr = wvlt_init_cond_trans(f0r,psirf);
fwci = wvlt_init_cond_trans(f0i,psiif);

hfunr = Kmat*hfltr;
hfuni = Kmat*hflti;
gfunr = Kmat*gfltr;
gfuni = Kmat*gflti;

hfltrr = real(Kmati*hfunr/KT);
hfltir = real(Kmati*hfuni/KT);
gfltrr = real(Kmati*gfunr/KT);
gfltir = real(Kmati*gfuni/KT);

dwtcr(lstp:rstp) = fwcr;
dwtci(lstp:rstp) = fwci;

dwtmode('per')

for jj=1:nlvls
    [ar,dr] = dwt(f0r,gfltrr,hfltrr);  
    [ai,di] = dwt(f0i,gfltir,hfltir);  
    
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
