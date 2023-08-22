function [SymbolsDetect] = deMod(IQin,Param)

switch Param.modulationFormat
    case {'dp-qpsk'}
        for Idx = 1:2
            IQin(:,Idx) = IQin(:,Idx)/rms(IQin(:,Idx));
        end
        SymbolsDetect = qamdemod(IQin,Param.M,'UnitAveragePower',true); % demodulation
    case {'dp-mpsk'}
        [SymbolsDetect, ~] = PSKdeMod(IQin,Param); % this is my function, can use Matlab insted
    case {'PAM'}
        SymbolsDetect = pamdemod(IQin,Param.M);
    otherwise
        error('Wrong modulation format specified');
end
end

