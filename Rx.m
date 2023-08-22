function [SigRx] = Rx(SigChannel,noise,Param)

if Param.OpticsEn
    %--- Photo-diode Operation
    SigRx = Param.PD.Resp*abs(SigChannel).^2;
    SigRx = SigRx - mean(SigRx); % reduce DC
else
    SigRx = SigChannel;
end
SigRx = SigRx + noise;         % Add noise


% Rx filter
SigRx = filter(Param.state.Brx,Param.state.Arx,SigRx); % Rx

%--- ADC filter
SigRx = filter(Param.state.Badc,Param.state.Aadc,SigRx); % Rx
% quantization noise is neglected 

end


