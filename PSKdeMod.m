function [symbolDet,IQdetect] = PSKdeMod(IQin,Param)
M = Param.M;
initialphase = Param.dataGenerator.initialphase;
IQin = IQin*exp(1i*(pi/M - initialphase)); % make the phase more convenient to partition 
symbolAmp = (0:M-1);
ThreshVecPhase = symbolAmp*2*pi/M ;
ThreshVecPhase(end+1) = 2*pi;

symbolDet = nan(size(IQin));
Angle =  mod(angle(IQin),2*pi); %[rad] range [0,2*pi]

for Idx = 1:M
    Index = (ThreshVecPhase(Idx) <= Angle) & (Angle < ThreshVecPhase(Idx+1));
    symbolDet(Index) = symbolAmp(Idx);
end
IQdetect = exp(1i*(2*pi*symbolDet/M + initialphase));

