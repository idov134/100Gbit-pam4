clc;
clear;
close all;
%..............................................................
% ---------Setting parameters
%..............................................................
Rb = 50e9;       % symbol rate
sps = 8;
M = 4;
NOS = 2^17;

[ Param ] = LoadDefaultParamN(NOS,Rb,sps,M);

% --- Add Local param
% CLLE
NumOfCTLE = 12;
Param.DCGain = 0:-1:-NumOfCTLE; %DC Gain. 
Param.PeakingGain = 0:1:NumOfCTLE; %Peaking Gain
Param.PeakingFrequency = Param.Rs/2; %Peaking Frequency

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

% --- Overide some state param
if 0
Param.state.Badc = 1; % disable ADC filtering at the Rx since the CTLE is added after RX filter but before the ADC filter
Param.state.Aadc = 1; % disable ADC filtering at the Rx since the CTLE is added after RX filter but before the ADC filter
end


%% --- start the Transmission
%............................................................
% --- Transmitter
%.........................................................

%--- Source generation and modulation-----------------
% Generate random data source to be transmitted
[ModulatedData, Param] = dataGenerator(Param);
symbols = Param.state.symbols;
SigTransmitted = TxPS(ModulatedData, Param);

%...............................................................
%-------- Channel and AWGN
%.............................................................
FIRlen = 50;


if 1
    [SigChannel,Param] = FilterS21(SigTransmitted, S21strc, FIRlen, Param); % 'SigChannel' is [1048576x1] 
    % symbols that were holded, passed through the Tx and the PCB.
    SigChannel = real(SigChannel); % For some reason it thinks that 'SigChannel' is complex (0.000000i).
else % We will enter this else if we don't have the S21 information. (but we have)
    load(['C:\Gilad\Edi\ImpuseReponse.mat'])
    Param.ChannelResponse = h;
    FsImp =  26.56e9*8;
    channel = interp1(Param.ChannelResponse,1:FsImp/Fs:length(Param.ChannelResponse),'cubic'); channel = channel/sum(channel);
    % channel = [1 ]; % for debuging
    %freqz(channel,1,1e4,Fs)
    SigChannel = filter(channel, 1, SigTransmitted);
end

% --- AWGN
rng(Param.Rx.NoiseSeed)
noise = randn(size(SigChannel))*sqrt(Fs/(2*F_cut_rx)); % random noise with the same length as 'SigChannel'
% *sqrt(Fs/(2*FpdCut)) - consider the Rx filter noise reduction.
sweepParam = SNR; % 0:5:35[dB].

BER = nan(length(sweepParam),1); % [8x1] Column vector that represents the BER for every SNR (initializing as NaN).


% ----- Penalty Calculations:
NumFfeTapsVec = [5, 4, 3]; % We will check the penalty for this values.
NumDfeTaps = 1;


for NumFfeTaps_i = NumFfeTapsVec 
    for s = 1:NumOfCTLE+1
        
        % --- sweep (the SNR)
        for Idx = 1:length(sweepParam) % Every iteration will give us the BER for different SNR.
        
            % --- the channel attenuation
            SigChannel = SigChannel/rms(SigChannel); % Normalizing the symbols [1048576x1] that were holded and passed through the Tx and the PCB.
            SigChannel = SigChannel*10.^(+SNR(Idx)/20); % The larger the SNR, the greater the attenuation.
        
            % --- receiver part
            SigRx = Rx(SigChannel,noise,Param); % 'SigRx' are the symbols that were holded passed through 
            % the Tx, PCB, Rx and ADC.
        
            % --- CTLE
%             if Idx == 0
             ImpusleResCTLE = AddCTLE_MinBER(SigRx,Param);
             CTLE = ImpusleResCTLE(:, s);
             SigRx = filter(CTLE,1,SigRx);

%             end

            if 0
                eyediagram(SigRx(103:1e3),sps)
            end
        
            % --- sample at Fs = Rs
            [SigRxUsmp,~] = Sampler(SigRx,sps);
            SigRxUsmp = SigRxUsmp/max(abs(SigRxUsmp))*max(ModulatedData); % normalizing
        
            % --- FFE/DFE
            [~,w,c,Error,dataOut] = DfeLms(SigRxUsmp,ModulatedData,NumFfeTaps_i,Mu,CursorPos,NumDfeTaps);
            %     load('w_and_c')
        
            % ---using the function with a fix FFE/DFE weights
            SigDfeOut = DfeLms(SigRxUsmp,dataOut,NumFfeTaps_i,Mu,CursorPos,NumDfeTaps,0,w,c,1);
% % % % % % % % % % % % % % %                 figure(55); hist(SigDfeOut,50); title('after DFE'); figure(56) ; hist(SigRxUsmp,50) ; title('before DFE');
        
            % --- deModulation
            DataRx = deMod(SigDfeOut,Param);
        
            %--- BER calculation
            [BER(Idx), ~] = getBER(DataRx,Param);
        
% % % % % % % % % % % % % % %                 figure(555);
% % % % % % % % % % % % % % %                 semilogy(SNR(1:Idx),BER(1:Idx));
        end
        All_BERS(:, NumFfeTaps_i, s) = BER; % For Penalty.
    end
end


% % % % % % % % % % % % % % % figure(555)
% % % % % % % % % % % % % % % ylabel('BER'); grid on
% % % % % % % % % % % % % % % xlabel('SNR [dB]')


% ---- Plot Penalty Results ----
BERPreferably = 1e-3;
for s = 1:NumOfCTLE+1
    for NumFfeTaps_i = NumFfeTapsVec
        [~, ~, ~, ~, penalty_specificBER_db(NumFfeTaps_i, s)] = FindPenalty(SNR, All_BERS(:, NumFfeTaps_i, s), SNR, All_BERS(:, NumFfeTapsVec(1), 1), BERPreferably);
    end
end
figure(1234)
plot(1:NumOfCTLE+1, penalty_specificBER_db(NumFfeTapsVec, 1:NumOfCTLE+1), '-*');
xlabel('CTLE configuration');
ylabel('Penalty [dB]');
title('CTLE - Penalty vs. CTLE configuration for BER = 1\times 10^{-3}');
legend('for 5 FFE Taps', 'for 4 FFE Taps', 'for 3 FFE Taps');
grid on;



% Plot best configurations of all the configurations of the CTLE:
for z = 1:NumOfCTLE+1
    figure(12345);
    N = length(ImpusleResCTLE(:, z));
    FreqResCTLE(:, z) = fftshift(fft(ImpusleResCTLE(:, z)));
    Freq = (0:N/2-1)*Fs/N;
    plot(Freq./1e9, db(abs(FreqResCTLE(1:N/2, z)))); % We plot all the CTLE's.
    hold on;
    title('All The Configurations of CTLE')
    xlabel('Frequency [GHz]');
    ylabel('Magnitude [dB]');
    legend('1st configuration', '2nd configuration', '3rd configuration', '4th configuration', '5th configuration',...
        '6th configuration', '7th configuration', '8th configuration', '10th configuration', '11th configuration',...
        '12th configuration', '13th configuration');
    grid on;
%     ylim([-60, 0]); % set the range of the y-axis in dB
%     yticks(-60:10:0); % set the tick labels for the y-axis
%     yticklabels(string(-60:10:0)+" dB"); % set the tick labels with "dB" units
end
hold off;



% Plot results
figure(1);
subplot(2,1,1)
semilogy((Error.^2));
title(['LMS Adaptation Learning Curve Using Mu = ', num2str(Mu)]);
xlabel('Iteration Number');
ylabel('Output Estimation Error in dB');
subplot(2,1,2)
stem(w); hold on ; stem(length(w)+1:length(w)+ length(c),c,'r') ; hold off
title('DFE taps')
% subplot(3,1,3)
% stem(real(channel)); title('Channel Taps')
