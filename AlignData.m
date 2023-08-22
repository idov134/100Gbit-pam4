function [Sig1_alinged,Sig2_alinged,idxOut] = AlignData(Sig1,Sig2, is_circ)
% align Sig1 with Sig2
if nargin<3
    is_circ=1;
end
xc=xcorr(Sig1,Sig2);
[~,idx] = max(xc);
idxOut = idx-length(Sig1);
if is_circ
    Sig1_alinged = circshift(Sig1, -idxOut);
    Sig2_alinged = Sig2;
else
    if idxOut > 0
        Sig1_alinged = Sig1(idxOut+1:end);
        Sig2_alinged = Sig2;
    else
        Sig2_alinged = Sig2(-idxOut+1:end);
        Sig1_alinged = Sig1;
    end
end
len=min(length(Sig1_alinged),length(Sig2_alinged));
Sig1_alinged = Sig1_alinged(1:len);
Sig2_alinged = Sig2_alinged(1:len);
end