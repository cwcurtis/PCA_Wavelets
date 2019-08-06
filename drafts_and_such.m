drafts

KTT = osamp*KT;
nones = (-1).^(-KTT/2+1:KTT);
gfltr = nones.*gfltr;
gflti = nones.*gflti;
hfltr = nones.*hfltr;
hflti = nones.*hflti;

omask = ones(KTT,1);
omask(2:2:KTT) = 0;
emask = 1 - omask;

gfltre = emask.*gfltr;
gfltro = omask.*gfltr;
gfltie = emask.*gflti;
gfltio = omask.*gflti;

hfltre = emask.*hfltr;
hfltro = omask.*hfltr;
hfltie = emask.*hflti;
hfltio = omask.*hflti;

hfunre = Kmat*hfltre;
hfunro = Kmat*hfltro;
hfunie = Kmat*hfltie;
hfunio = Kmat*hfltio;

gfunre = Kmat*gfltre;
gfunro = Kmat*gfltro;
gfunie = Kmat*gfltie;
gfunio = Kmat*gfltio;

ffunre = Kmat*(emask.*f0r);
ffunro = Kmat*(omask.*f0r);
ffunie = Kmat*(emask.*f0i);
ffunio = Kmat*(omask.*f0i);
