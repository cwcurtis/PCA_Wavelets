% This function generates the 2pi - periodization of the input function
% f(x).

function ffold = gam_comp(dx,fin,gin,Ncp,Xmesh)
    ffold = zeros(length(fin),1);
    for jj = -Ncp:Ncp
       ffold = ffold + dx*exp(-2*pi*1i*jj*Xmesh).*fft(fin);
    end
end