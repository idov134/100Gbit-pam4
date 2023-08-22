function [F_rx_BW] = GetRxBW(Param)
% Get the Rx inclusive BW (3dB point)
b = conv(Param.state.Badc,Param.state.Brx);
a = conv(Param.state.Aadc,Param.state.Arx);
[h,F] = freqz(b,a,1e4,Param.Fs); % freq response of the adc and the rx 
% plot(F,db(h))
[~,IndexCut] = min(abs(db(h) + 3));
F_rx_BW = F(IndexCut);

end