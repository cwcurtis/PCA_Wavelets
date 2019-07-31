function rhs = nonlinearity(K,eta,q,G0,ep,sig,Kmesh)

KT = 2*K;
% Find the wave numbers to implement the 2/3 de-aliasing throughout
Kc = floor(2*K/3);
Kuc = KT-Kc+1;
Kc = Kc+1;
Dx = 1i*Kmesh;
Dx2 = -Kmesh.^2;
G0op = abs(Dx);
Hop = 1i*sign(Kmesh);

eta(Kc:Kuc)=0;
q(Kc:Kuc)=0;
etax = real(ifft(Dx.*eta)); 
etaxx = real(ifft(Dx2.*eta)); 
qx = real(ifft(Dx.*q)); 
etap = real(ifft(eta));

G1 = real(ifft( -G0op.*fft(etap.*G0) - Dx.*fft(etap.*qx) ));
G2 = real(ifft( -G0op.*fft(etap.*G1) + .5*Dx2.*( fft(etap.^2.*G0) + Hop.*fft(etap.^2.*qx) ) ));

rhs1 = fft( ep*G1 + ep^2*G2 );
rhs2 = fft( .5*ep*(-qx.^2 + G0.^2 + 2*ep*G0.*G1) + ep^2*(etax.*qx.*G0-3/2*sig*etax.^2.*etaxx) );

rhs1(Kc:Kuc) = 0;
rhs2(Kc:Kuc) = 0;

rhs = [rhs1;rhs2];
