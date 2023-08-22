clc;
clear;
close all;
%..............................................................
% ---------Setting parameters
%..............................................................
Rb = 100e9;       % symbol rate
sps = 8;
M = 4;
NOS = 2^17;

[ Param ] = LoadDefaultParamN(NOS,Rb,sps,M);

% Changing the number of taps of the DFE/FFE to see what will heppend:
Param.DFE.CursorPos = [4, 4, 4, 4];
Param.DFE.numDfeTaps = [3, 3, 0, 3];
Param.DFE.numFfeTaps = [10, 4, 10, 3];


%---------------------Setting other Parameters------------------

%--- General

Param.state.frame_inx = 1; %future purpose

%--- Channel
SNR = Param.SNR;

% --- FFE
% Setting Other paramters
CursorPos = Param.DFE.CursorPos;
numDfeTaps = Param.DFE.numDfeTaps;
numFfeTaps = Param.DFE.numFfeTaps;	 % channel order
Mu = Param.DFE.Mu;          % iteration step size
% Mu = 0.05;

%--- overriding paramters
Param.OpticsEn = false;

%--- Building  other Parameters
Fs = Param.Fs; % sample rate, time two since the spectrum is doubled (to make the ifft real)
% --- handle PCB
Param.CTLEen = false; % enable/dis CTLE
filename = 'THRU_LL_Rx1_L_Diff.mat';
S21strc = load(filename);
%rfplot(S21strc)

Param = GetFilterCoef(Param); %   Getiing the filters Coefs (b and a)
F_cut_rx = GetRxBW(Param); % Get the Rx inclusive BW


%% --- start the Transmission
%............................................................
% --- Transmitter
%.........................................................

%--- Source generation and modulation-----------------
% Generate random data source to be transmitted
[ModulatedData, Param] = dataGenerator(Param);
symbols = Param.state.symbols;
SigTransmitted = TxPS(ModulatedData, Param);


% --- eyediagram for the transmitted symbols (after the ZOH and Tx):
eyediagram(SigTransmitted(103:1e3),sps)
title('Transmitted symbols (after the ZOH and Tx)')


%...............................................................
%-------- Channel and AWGN
%.............................................................
FIRlen = 50;

[SigChannel,Param] = FilterS21(SigTransmitted, S21strc, FIRlen, Param); % 'SigChannel' is [1048576x1] 
% symbols that were holded, passed through the Tx and the PCB.
% We will enter this 'else' if we don't have the S21 information. (but we have)

% --- eyediagram for the symbols after the pcb (without noise) (after the ZOH, Tx and pcb):
eyediagram(SigChannel(103:1e3),sps)
title('symbols after the pcb (without noise)')

% --- AWGN
rng(Param.Rx.NoiseSeed)
noise = randn(size(SigChannel))*sqrt(Fs/(2*F_cut_rx)); % random noise with the same length as 'SigChannel'
% *sqrt(Fs/(2*FpdCut)) - consider the Rx filter noise reduction.
sweepParam = SNR; % 0:5:35[dB].

BER = nan(length(sweepParam),1); % [8x1] Column vector that represents the BER for every SNR (initializing as NaN).

% --- sweep (the SNR)
A = 100; % For the figures.
for i = 1:4
    for Idx = 1:length(sweepParam) % Every iteration will give us the BER for different SNR.
    
        % --- the channel attenuation
        SigChannel = SigChannel/rms(SigChannel); % Normalizing the symbols [1048576x1] that were holded and passed through the Tx and the PCB.
        SigChannel = SigChannel*10.^(+SNR(Idx)/20); % The larger the SNR, the greater the attenuation.

        % --- receiver part
        SigRx = Rx(SigChannel,noise,Param); % 'SigRx' are the symbols that were holded passed through 
        % the Tx, PCB, Rx and ADC.


        % --- eyediagram for the symbols after all the system (after the ZOH, Tx, PCB, Rx and ADC):
        eyediagram(SigRx(103:1e3),sps)
        title(['symbols after all the system SNR = ', num2str(SNR(Idx))])

        % --- sample at Fs = Rs
        [SigRxUsmp,~] = Sampler(SigRx,sps);
        SigRxUsmp = SigRxUsmp/max(abs(SigRxUsmp))*max(ModulatedData); % normalizing

        % --- FFE/DFE
        [~,w,c,Error,dataOut] = DfeLms(SigRxUsmp,ModulatedData,numFfeTaps(i),Mu,CursorPos(i),numDfeTaps(i));
        %     load('w_and_c')

        % ---using the function with a fix FFE/DFE weights
        SigDfeOut = DfeLms(SigRxUsmp,dataOut,numFfeTaps(i),Mu,CursorPos(i),numDfeTaps(i),0,w,c,1);
        figure(Idx*10); hist(SigDfeOut,50); title(['after DFE SNR = ', num2str(SNR(Idx))]); % The optimum sample (symbols) after the DFE FFE. 
        figure(A) ; hist(SigRxUsmp,50) ; title(['before DFE SNR = ', num2str(SNR(Idx))]); % The optimum sample before (symbols) the DFE FFE.
        A = A + 100;
    

        % --- deModulation
        DataRx = deMod(SigDfeOut,Param);

        %--- BER calculation
        [BER(Idx), ~] = getBER(DataRx,Param);

        figure(555);
        semilogy(SNR(1:Idx),BER(1:Idx));
    end
    hold on;
end
hold off;
figure(555)
ylabel('BER'); grid on
xlabel('SNR [dB]')
title('BER')
legend('FFE taps = 10 & DFE taps = 3', 'FFE taps = 4 & DFE taps = 3', 'FFE taps = 10 & DFE taps = 0', 'FFE taps = 3 & DFE taps = 3');

% Plot results

figure(1000);
% semilogy((Error(11:end).^2));
plot((Error(11:end)));
envelope((Error(numFfeTaps+1:end)), 30, 'peak');
xlim([0, 1e4]);
ylim([-4, 4]);
title(['LMS Adaptation Learning Curve Using Mu = ', num2str(Mu)]);
xlabel('Iteration Number');
ylabel('Output Estimation Error in dB');
figure(1100);
subplot(2, 1, 1)
stem(w); hold on ; stem(length(w)+1:length(w)+ length(c),c,'r') ;
legend('FFE', 'DFE')
hold off
title('DFE and FFE Taps')
subplot(2, 1, 2)
stem(real(Param.state.PCB.FIR)); title('Channel Taps')
xlim([0, 10]);
% We don't sure if
% stem of the zeros of the pcb will give us the taps of the
% pcb!!!!!!!!!!!!!!!


figure(2000)
[H,F] = freqz(Param.state.PCB.FIR, 1, 1e4, Param.state.PCB.Fs); % The frequency response of the channel.
plot(F/1e9,db(H))
title('Frequency response of the channel');
xlabel('Frequency [GHz]');
ylabel('Magnitude [dB]');
grid

% % figure(2100)
% % rfplot(S21strc.SparamDiff) % The frequency response of the channel (all the S param of the channel).
% % title('Frequency response of the channel');
