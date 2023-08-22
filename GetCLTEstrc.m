function CTLEstrc = GetCLTEstrc(Param)

% This function gets:
% *Param - The struct with all our parameters.

% And returns:
% * 'CTLEstrc' - an CTLE object that contain all the parameters that were
% defined.

SymbolTime = 1/Param.Rs;
DCGain = 0:-1:-26; % Gain at zero frequency for the CTLE transfer function. (!same length as 'PeakingGain'!).
PeakingGain = 0:26; % The difference between the maximum gain and the gain in the beginning (0 Hz).
PeakingFrequency = Param.Rs*0.7; % The peak gain will be at 0.7*Rs (=0.7*25GHz).
% Will be converted to match the length of 'DCGain' and 'PeakingGain' by scalar expansion.
dt = 1/Param.Fs; % (=1/200G).

CTLEstrc = serdes.CTLE('SymbolTime',SymbolTime,'SampleInterval',dt,...
    'Mode',2,'WaveType','Impulse',...
    'DCGain',DCGain,'PeakingGain',PeakingGain,...
    'PeakingFrequency',PeakingFrequency); 
% *An impulse response input signal will be.
% *CTLE Mode is 'adapt' (=2) - The serdes.CTLE determines the CTLE transfer 
% function to maximize the performance metric as specified by the PerformanceCriteria 
% property and applies the transfer function to the input waveform for time domain simulation.
end

