function ImpusleResCTLE = AddCTLE(SigIn,Param)

% The function gets:
% *SigIn - The symbols that were holded, passed through the Tx, PCB, Rx
% and ADC (='SigRx').
% *Param - Our struct with all the parameters.

% This function returns:
% *ImpusleResCTLE - the best impulse response of the CTLE from all the
% configurations of CTLE that were tested.


SymbolTime = 1/Param.Rs; % =4e-11.
sps = Param.sps; % =8
dT = 1/Param.Fs; % =1/200G = 5e-12.

CTLE1 = serdes.CTLE(...
    'SymbolTime',SymbolTime,...     %Duration of a single symbol 
    'SampleInterval',dT,...         %time step size
    'DCGain',Param.DCGain,...       %DC Gain (Gain at zero frequency !same length as 'PeakingGain'!) (=0:-1:-12).
    'PeakingGain',Param.PeakingGain,...  %Peaking Gain (The difference between the maximum gain and the gain in the beginning) (=0:1:12).
    'PeakingFrequency',Param.PeakingFrequency,...  %Peaking Frequency (The peak gain will be at..) (=1.25e10) (will become a vector as 'PeakingGain').
    'Mode',1);                  %Mode is fixed
% For each configuration of the CTLE, stimulate the CTLE with an ideal step response excitation to extract the reference CTLE step responses and observe the output waveforms.
% The defult of 'ConfigSelect' is 0 (if =5 than the 6's configuration will be selected).

stimulus = [1 ; zeros(10*sps,1)]; % impulse
% stimulus(1:SamplesPerSymbol) = 0;  (column vector of '1' and 80 '0' (81x1)).

numberOfConfig = CTLE1.ConfigCount; % number of configurations of CTLE (=13).

N = min(2e3,length(SigIn)); % N = min(2e3, 1048576) = 2000.
stepResponse = zeros(length(stimulus),numberOfConfig); % [81x13] matrix of zeros (every column will represent the impulse response of a specific configuration of CTLE (the input is impulse response)).
Peak2Avg = nan(numberOfConfig,1); % 13x1 vector of NaN.
SigOutCTLE = nan(N,numberOfConfig); % [2000x13] matrix of NaN. (every column represents the output of a specific impulse reponse of CTLE when the input is 'SigIn(1:N)' (holded symbols after Tx, PCB, Rx)).
for Idx = 1:numberOfConfig % Every iteration is different configuration.
    CTLE1.ConfigSelect = Idx-1; % The configuration that was selected.
    release(CTLE1); % Allows to change properties and input characteristics in the object 'CTLE1'.
    stepResponse(:,Idx) = CTLE1(stimulus); % Every column represents a impulse response of different configuration of CTLE.
    SigOutCTLE(:,Idx) = filter(stepResponse(:,Idx),1,SigIn(1:N)); % conv with all CTLE options - 
    % (every column represents the output of a specific impulse reponse of CTLE when the input is 'SigIn(1:N)' (holded symbols after Tx, PCB, Rx)).
    % = The symbols passed through different configurations of CTLE.
    Peak2Avg(Idx) = db(max(SigOutCTLE(:,Idx)/rms(SigOutCTLE(:,Idx)))); % The largest normalized symbol (in dB).
end


% --- Choose the best CTLE
% The best CTLE will be the CTLE that will give the smallest (between all
% the configuration) largest normalized symbol (in dB).
[~,Idxmin] = min(Peak2Avg); % The 4'th configuration of CTLE is the best.
ImpusleResCTLE = stepResponse(:,Idxmin); % The impulse response of the 4'th CTLE. (81x1).

if 0
eyediagram(filter(ImpusleResCTLE,1,SigIn(1:N)),sps)
figure,
t1 = dT*(0:size(stepResponse,1)-1);
plot(t1,stepResponse)
xlim([0, 1e-10]);
xlabel('time [s]'),ylabel('[V]')
title('CTLE Step responses')
legend(cellstr(num2str((0:(numberOfConfig-1))')))
grid on
end

% C:\Users\Gilad\Documents\MATLAB\Examples\R2021b\serdes\ConvertSParameterToImpulseResponseForSerDesToolboxExample
