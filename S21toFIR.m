function [FIR] = S21toFIR(S21, F, FirLen)
%input  S21 - freq. response in to out
% F - frequency coresponding to S21
% FilLen - the FIR length w/o zero padding - the resu;

% This function gets:
% *S21 - The 4000x1 vector that represents the s21 parameters (the 2nd row, the 1st column and all
% the depthes from the tensor that was loaded).
% *F - 4000x1 vector that represents the frequency (from 0.1G to 40G).
% *FirLen - (='FIRlen'=50).

% and returns:
% *FIR - The ifft of the interpolated (to be 25x1 (from 4001x1)) S21 (='S21new') and the double side
% (perfect symmetry) of 'S21new' (that the first angle is 0), and than
% normalized. (49x1 vector).
% @@@@@@@@@@@@@@IN EASY TERMS - 'FIR' - is the zeros of the digital filter that
% represents the PCB response (freqz(FIR, 1) will give us the frequency
% response)@@@@@@@@@@@@@@@@@@
 
if F(1) ~= 0 % (F(1) = 1e7).
    tmp = unwrap(angle(S21(1:2))); % The phases of the first 2 elements in S21 (the jump will be no more than pi).
    Phase = tmp(1) - diff(tmp); % phase estimation S21(F=0). (we assume that the diff between the 0 to the 1st = diff between the 1st and the 2nd). 
    S21 = [abs(S21(1))*exp(1i*Phase); S21]; % Adding the estimated S21 element to be the first in the S21 vector (now it is 4001x1).
    F = [0;F]; % The frequency of the first S21 element is 0 (now it is 4001x1).
end 
FirLen = FirLen/2; % (=25)
if FirLen < length(S21) % 25 < 4001.
    q = length(S21)/FirLen; % (=4001/25)
else
    sprintf('Fir Length was not incresed to %d',FirLen)
end
S21new = interp1((S21),(1:q:length(F))','spline'); % 'S21new' will have the same behavior as 'S21' and its' 'x' will be "(1:q:length(F))'". (25x1 complex vector).
S21new = S21new.*exp(-1i*angle(S21new(1))); % to set the phase of S21new(1) to 0 (extracting 'angle(S21new(1))' from all the elements)(but it changes also a bit of the other values).
S21_doubleSide = [S21new; flipud(conj(S21new(2:end)))]; % Creating the double side (perfect symmetry) from 25:49 (flip up to down on the conj). (49x1 complex vector).                
FIR = (ifft(S21_doubleSide));
FIR = FIR/sum(FIR);

end