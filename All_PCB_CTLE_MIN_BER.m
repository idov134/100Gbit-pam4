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
% % % filename = 'THRU_LL_Rx1_L_Diff.mat';
% % % S21strc = load(filename);
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

% ----- For All The PCBs:
filenames = ["31dB_at_26p56G_stripline.mat", ...
    "26dB_at_26p56G_stripline.mat", "21dB_at_26p56G_stripline.mat", ...
    "15p7dB_at_26p56G_stripline.mat", "10p4dB_at_26p56G_stripline.mat", "THRU_LL_Rx1_L_Diff.mat"];

% % % % filenames = ["31dB_at_26p56G_stripline.mat", ...
% % % %     "26dB_at_26p56G_stripline.mat"];


s = 1; % Parameter for the CTLE plots.
v = 1; % Parameter for the All BERs.
for names = filenames
    S21strc = load(names);
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
    NumFfeTapsRef = Param.DFE.numFfeTaps; % =16.
    NumDfeTapsRef = Param.DFE.numDfeTaps; % =4.
    
%     NumFfeTapsVec = [NumFfeTapsRef, 4, 7, 10, 13, 16]; % We will check the penalty for this values.
%     NumDfeTapsVec = [NumDfeTapsRef, 3, 4]; % We will check the penalty for this values.

    NumFfeTapsVec = [NumFfeTapsRef, 13, 10, 7, 4, 3]; % We will check the penalty for this values.
    NumDfeTapsVec = [NumDfeTapsRef, 3, 2]; % We will check the penalty for this values.
    
    for NumDfeTaps_j = NumDfeTapsVec
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
        
%                     if NumDfeTaps_j == 3 && NumFfeTaps_i == 7 % See CTLE best conf for 3 DFETaps and 7 FFETaps.
%                         BestCTLEImpulseRes(:, s) = ImpusleResCTLE(:, BestCTLE37(v)); % 81x48 matrix - the first 8 columns is for the first PCB, and each column is another SNR (same for the second 8...)
%                         N = length(CTLE);
%                         FreqResCTLE(:, s) = fftshift(fft(BestCTLEImpulseRes(:, s)));
%                         Freq = (0:N/2-1)*Fs/N;
%                     end
        
                    if 0
                        eyediagram(SigRx(103:1e3),sps)
                    end
                
                    % --- sample at Fs = Rs
                    [SigRxUsmp,~] = Sampler(SigRx,sps);
                    SigRxUsmp = SigRxUsmp/max(abs(SigRxUsmp))*max(ModulatedData); % normalizing
                
                    % --- FFE/DFE
                    [~,w,c,Error,dataOut] = DfeLms(SigRxUsmp,ModulatedData,NumFfeTaps_i,Mu,CursorPos,NumDfeTaps_j);
                    %     load('w_and_c')
                
                    % ---using the function with a fix FFE/DFE weights
                    SigDfeOut = DfeLms(SigRxUsmp,dataOut,NumFfeTaps_i,Mu,CursorPos,NumDfeTaps_j,0,w,c,1);
    % % % % % % % % % % % % % % % % % %                 figure(55); hist(SigDfeOut,50); title('after DFE'); figure(56) ; hist(SigRxUsmp,50) ; title('before DFE');
                
                    % --- deModulation
                    DataRx = deMod(SigDfeOut,Param);
                
                    %--- BER calculation
                    [BER(Idx), ~] = getBER(DataRx,Param);
                
    % % % % % % % % % % % % % % % % % % % % %                 figure(555);
    % % % % % % % % % % % % % % % % % % % % %                 semilogy(SNR(1:Idx),BER(1:Idx));
                end
                CTLEBERs(:, s, v) = BER;
            end
            BERsSNR20(:, v) = CTLEBERs(5, :, v); % BER matrix for SNR = 20 for every config of CTLE (every column is different PCB and every row is different conf).

            % --- Finding The CTLE That Gives Min BER (for every PCB):
            [~, BestCTLEIdx(v)] = min(BERsSNR20(:, v));

            if NumDfeTaps_j == 3 && NumFfeTaps_i == 7 % See CTLE best conf for 3 DFETaps and 7 FFETaps.
                BestCTLE37(v) = BestCTLEIdx(v);
            end

            All_BERS(:, NumFfeTaps_i, NumDfeTaps_j, v) = CTLEBERs(:, BestCTLEIdx(v), v); % A 4 dimensional tensor the rows is the individual BER for specific SNR, the column in for specific number of FFE taps, the depthes is for specific DFE taps, and the 4's dimension is for specific PCB.
        end
    end
    v = v + 1; % In order to change PCB.
end

% % % % % % % % % % % % % % % % % % % % figure(555)
% % % % % % % % % % % % % % % % % % % % ylabel('BER'); grid on
% % % % % % % % % % % % % % % % % % % % xlabel('SNR [dB]')


% ---- Plot BER as Function of Configuration Number (of CTLE) for 2 DFE taps and 3 FFE taps ----
for v = 1:6
    figure(1111*v)
    stem(1:NumOfCTLE+1, BERsSNR20(:, v));
    title('The BER for every CTLE at SNR = 20');
    xlabel('Number of Configuration');
    ylabel('BER');
end

% ---Plot the BER as function of PCB (for the best CTLE) at SNR = 20 for 2 DFE taps and 3 FFE taps--- 
figure(999)
stem(1:6, min(BERsSNR20));
title('The BER for every PCB at SNR = 20 dB');
xlabel('PCB');
xlim([0, 7]);
ylabel('BER');

% ---- Plot Penalty Results ----
BERPreferably = 1e-3;
for v = 1:6
    for NumDfeTaps_j = NumDfeTapsVec
        for NumFfeTaps_i = NumFfeTapsVec
            [~, ~, ~, ~, penalty_specificBER_db(NumFfeTaps_i, NumDfeTaps_j, v)] = FindPenalty(SNR, All_BERS(:, NumFfeTaps_i, NumDfeTaps_j, v), SNR, All_BERS(:, NumFfeTapsVec(1), NumDfeTapsVec(1), v), BERPreferably);  
        end
    end
    figure(111*v)
    plot(NumFfeTapsVec(1:end), penalty_specificBER_db(NumFfeTapsVec, NumDfeTapsVec, v), '-*');
    set(gca, 'XDir', 'reverse');
    hold on;
    xlabel('Number of FFE Taps');
    ylabel('Penalty [dB]');
    title('FFE - Penalty vs. Number of Taps for BER = 1\times 10^{-3}');
    legend('for 4 DFE Taps', 'for 3 DFE Taps', 'for 2 DFE Taps');
    grid on;
end
hold off;



% --- Plot the penalty as a function of the PCB for 16 taps of FFE and 4 taps of DFE (the ref is the 5'th (shortest) PCB)---
figure(666);
for v = 1:6
    [~, ~, ~, ~, penalty_pcb(v)] = FindPenalty(SNR, All_BERS(:, 16, 4, v), SNR, All_BERS(:, 16, 4, 6), 1e-2);
end
plot(1:6, penalty_pcb(1:6), '--o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerEdgeColor','b', 'Color', 'k');
set(gca, 'XDir', 'reverse');
hold on;
xlabel(' The PCB');
ylabel('Penalty [dB]');
title('Penalty vs. PCB for BER = 1\times 10^{-2}');
grid on;




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


% Plot best configurations of CTLE for every PCB (for 7 FFE taps and 3 DFE taps):
for z = 1:6
    figure(12345);
    BestCTLEImpulseRes(:, z) = ImpusleResCTLE(:, BestCTLE37(z));
    N = length(BestCTLEImpulseRes(:, z));
    FreqResCTLE(:, z) = fftshift(fft(BestCTLEImpulseRes(:, z)));
    Freq = (0:N/2-1)*Fs/N;
    plot(Freq./1e9, db(abs(FreqResCTLE(1:N/2, z)))); % We choose the CTLE for 7 FFE taps and 3 DFE taps AND its' BER is closest to 1e-3.
    hold on;
    title('The CTLE for 7 FFE taps and 3 DFE taps')
    xlabel('Frequency [GHz]');
    ylabel('Magnitude [dB]');
    legend('Best CTLE for 31dB-at-26p56G-stripline', 'Best CTLE for 26dB-at-26p56G-stripline', ...
        'Best CTLE for 21dB-at-26p56G-stripline', 'Best CTLE for 15p7dB-at-26p56G-stripline',...
        'Best CTLE for 10p4dB-at-26p56G-stripline', 'Best CTLE for THRU-LL-Rx1-L-Diff');
    grid on;
%     ylim([-60, 0]); % set the range of the y-axis in dB
%     yticks(-60:10:0); % set the tick labels for the y-axis
%     yticklabels(string(-60:10:0)+" dB"); % set the tick labels with "dB" units
end
hold off;

