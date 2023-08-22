function [ Param ] = LoadDefaultParamN(NOS,Rb,sps,M)

% Constant system % Parameters--------------------------------------------------------------

Param.NOS = NOS; % number of bits to process
Param.M = M;
Param.sps = sps; % oversampling factor:
Param.Rb = Rb;
Param.Rs = Param.Rb/log2(M);
Fs = Param.Rs*sps;
Param.Fs = Fs;
Param.NOB = NOS*log2(M);
Param.NumOfTx = 1; % number of Txs (for coherent set 2) 
if (round(log2(NOS)) - log2(NOS)) ~= 0
    error('NOS should be power of two')
end

%--- General
Param.OpticsEn = true; % enable optical channel
Param.FilterDis = false; % DAC filter Disable
Param.modulationFormat = 'PAM'; % 'dp-qpsk' 'dp-mpsk' 'PAM'


% --- seeds
SeedBasic = 200;
Param.dataGen.dataSeed = SeedBasic;
Param.Rx.PD.ThermalNoiseSeed = [SeedBasic + 100];%, SeedBasic + 200;SeedBasic + 300, SeedBasic + 400]; % per polarization and per detector side 

% data generator
Param.dataGenerator.PRBS_enable = false;
Param.dataGen.dataSeed

% --- Tx
Param.Tx.ShaperType = 'ZOH';     % option 'RC' ir 'ZOH' Raised cosine (RC) or zero order hold

% --- Modulator (used if Param.Tx.ShaperType == 'ZOH')
Param.Tx.FilterType = 'Bessel'; % Bessel or Butter
Param.Tx.FilterOrder = 2; 
Param.Tx.Fcut = Param.Rs*0.8;

%--- Channel
Param.ChannelResponse =  [ 0.05 1 0.02 ]; % this is an arbitrary choice
% Param.x_channel  = 0:1/Param.Rs:length(Param.ChannelResponse)/Param.Rs - 1/Param.Rs;

Param.OIP = (-25:1:-17)+8;       %[dBm] optical input power (optical scenarion
Param.SNR = 0:5:35;             % SNR Electrical scenarion 

% --- Rx/Photodiode
Param.PD.Irn = 10e-12;  %[A/sqrt(Hz)] input referred noise 
Param.PD.Resp = 0.95;    %[A/W] Photodiode Responsivity

Param.Rx.FilterOrder = 4; 
Param.Rx.Fcut = Param.Rs*0.7;  %[Hz] Photodiode filter frequency cutoff
Param.Rx.NoiseSeed = 30;


% --- ADC
Param.ADC.FilterOrder = 4;    % ADC filter order
Param.ADC.Fcut = Param.Rs*0.7; %[Hz] ADC filter frequency cutoff


%% Not in use yet 
% --- DAC
Param.DAC.Fcut = Param.Rs*0.7; %[Hz] DAC filter frequency cutoff
Param.DAC.FilterOrder = 4;    % DAC filter order

Param.Dac.FilterType = 'Bessel';

% --- Photodiode 
Param.PD.FilterType = 'Butter'; % Bessel; Butter


%--- DFE parameters
Param.DFE.CursorPos = 4;
Param.DFE.numDfeTaps = 4;
Param.DFE.numFfeTaps = 16;	
Param.DFE.Mu = 0.005;        % iteration step size

%--- optical Channel
FiberLength = 2; %[km]
Param.Fiber.CD = FiberLength*17e-3;      % [sec/m] D = 17ps/nm/km (=17e-3 sec/m/km) at 1550 nm  
Param.Fiber.DGD  = 0e-12;    % seconds; PMD differentioal group delay [s]
Param.Fiber.Theta = pi/4;     % the angle between the reference polarization and the PSP of the Fiber
Param.Fiber.Lambda = 1550*1e-9;  % meters


%---
end

