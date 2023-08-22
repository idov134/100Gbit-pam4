
% 2nd order bessel filter
function [bd,ad] = bessel(n,Wn)
%Wn is the 3dB point in freqz plot (like butter)

%3dB point edjustment (empiric numbers)
factor_vec=[1 1.248 1.363 1.457 1.545 1.63 1.71 1.785];
factor=factor_vec(n);
Wn=Wn*factor;


[z,p,k] = besselap(n);
% [z,p,k] = buttap(n);
[b,a] = zp2tf(z,p,k);
arg=1/(2*tan(Wn*pi/2));
[bd,ad] = bilinear (b,a,arg);

if nargout==0
    freqz(bd,ad,1000);
end

end




% if n==1;
%     factor=1;
% elseif n==2;
%     factor=1.248;
% elseif n==3;
%     factor=1.363;
% elseif n==4
%     factor=1.457;
% elseif n==5
%     factor=1.545;
% elseif n==6
%     factor=1.63;
% elseif n==7
%     factor=1.71;
% elseif n==8
%     factor=1.785;
% end


