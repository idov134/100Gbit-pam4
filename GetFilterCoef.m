function Param = GetFilterCoef(Param)
% -- Get Filter coefs
Fs = Param.Fs;
[Param.state.Badc, Param.state.Aadc] = bessel(Param.ADC.FilterOrder, Param.ADC.Fcut/Fs*2);
[Param.state.Brx,Param.state.Arx] = bessel(Param.Rx.FilterOrder, Param.Rx.Fcut/Fs*2);

if strcmp(Param.Tx.ShaperType,'ZOH')
    switch Param.Tx.FilterType
        case 'Bessel'
            [Param.state.Btx,Param.state.Atx] = bessel(Param.Tx.FilterOrder,Param.Tx.Fcut/Param.Fs*2);
        case 'Butter'
            [Param.state.Btx,Param.state.Atx] = butter(Param.Tx.FilterOrder,Param.Tx.Fcut/Param.Fs*2);
    end
end
if Param.FilterDis == true ; Param.state.Badc = 1; Param.state.Aadc = 1; Param.state.Bpd = 1; Param.state.Apd = 1; end %Bpd = 1; Apd = 1;end

    
end
