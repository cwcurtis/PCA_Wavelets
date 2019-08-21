function [rphi,rphif] = phirscl(vphi,dx)
    
    fphi = fft(vphi);
    rphif = sqrt(dx)*exp(1i*angle(fphi));
    rphi = 1/dx*ifft(rphif);
    