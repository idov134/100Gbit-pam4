function [IQout, TxState] = TxPS(IQin,Param)

% Over sampling and TX filter pulse shaping
% Gilad

sps = Param.sps;
beta = 0.2; % should be used in loadparam
span = 10;  % should be used in loadparam

IQosp = upsample(IQin,sps);
switch Param.Tx.ShaperType
    case 'RC' % Raised Cosine
        b = rcosdesign(beta,span,sps);
    case 'ZOH' % Zero Order Hold
        b = ones(1,sps);
end

IQout = filter(b,1,IQosp);

% Adding TX filter (DAC/modulator amplifier e.t.c)  only for the ZOH case!
if strcmp(Param.Tx.ShaperType,'ZOH')
    IQout = filter(Param.state.Btx,Param.state.Atx,IQout);
end

TxState = []; % temporary
if 0
    eyediagram(IQout(200:1e3),2*sps)
end
end