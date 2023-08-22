function [IQout, Param] = dataGenerator(Param)
% Data generator 
% Inputs:
% sim structure 
% output:
% data generator state structure
% IQ - symbol in IQ format (complex) at one sample per symbol 

NOB = Param.NOB;
M = Param.M;
NumOfTx = Param.NumOfTx; 

if Param.dataGenerator.PRBS_enable
    
    bits = PRBS_generator(NOB,Param.data_generator.PRBS_order,lane_count);
    % Convert bits to symbols
    %          bits(1:13,1)=[0 0 0 0 0 0 1 0 0 0 0 0 0];
    symbols = bin2dec(num2str(bits));
else
    %Generate symbols {0,1,...,M-1}
    if(Param.state.frame_inx == 1)
        Param.state.DataState = rng(Param.dataGen.dataSeed,'v5uniform');
    end
    rng(Param.state.DataState)
    symbols = randi([0,M-1],Param.NOS,NumOfTx);
    Param.state.DataState = rng;
    
    symbols(1,1:NumOfTx) = 0; 
    symbols(end,1:NumOfTx) = 0;

    %Convert symbols to bits
%     bits = dec2bin(symbols(:)) - 48;
end
% dataGeneratorState.bits = bits;
% dataGeneratorState.symbols = symbols;

% Generate complex symbols
switch Param.modulationFormat
    case {'dp-qpsk'}
        IQout = qammod(symbols,M); % modulation
    case {'dp-mpsk'}
        IQout = exp(1i*(2*pi*symbols/M + Param.dataGenerator.initialphase));
    case {'PAM'}
        IQout = real(pammod(symbols,M));
    otherwise
        error('Wrong modulation format specified!');
end
Param.state.symbols = symbols;
if 0
    plot(real(IQout(10:end,1)),imag(IQout(10:end,1)),'o');shg
    xlabel('In-Phase') ; ylabel('Quadrature') ; grid
    title(sprintf('Data Constellation M = %d',M))
end



