clc;
clear;
close all;
%..............................................................
% ---------Setting parameters
%..............................................................
Rb = 50e9;       % Bit rate
sps = 8;
M = 4;
NOS = 2^17;

[ Param ] = LoadDefaultParamN(NOS,Rb,sps,M);

% --- Add Local param
% CTLE
NumOfCTLE = 12; % Number of Taps
Param.DCGain = 0:-1:-NumOfCTLE; % DC Gain
Param.PeakingGain = 0:NumOfCTLE; % Peaking Gain
Param.PeakingFrequency = Param.Rs/2; % Peaking Frequency

%---------------------Setting other Parameters------------------

%--- General

Param.state.frame_inx = 1; % future purpose

%--- Channel
SNR = Param.SNR;

% --- FFE
% Setting Other paramters
CursorPos = Param.DFE.CursorPos;
%numDfeTaps = Param.DFE.numDfeTaps;
%numFfeTaps = Param.DFE.numFfeTaps;	 
Mu = Param.DFE.Mu;                   % iteration step size

Checker = input('Choose FFE penalty - 1 or DFE penalty - 2 or Regular use - 3: '); 
switch Checker 
    case 1 % FFE penalty
        numFfeTaps_list = [8,10,12,14,16];
        numDfeTaps_list = [2,3,4,5,6];
    case 2 % DFE penalty
        numDfeTaps_list = [2,3,4,5,6];
        numFfeTaps_list = [4,6,7,8,10,12,14,16];
    case 3 % Regular use 
%         numFfeTaps_list = [10,16,9,11,8,7];
%         numDfeTaps_list = [1,2,5,5,6,6];
          numFfeTaps_list = 1:30;
          numDfeTaps_list = 1:30;

end

%--- overriding paramters
Param.OpticsEn = false;

%--- Building  other Parameters
Fs = Param.Fs; % sample rate, time two since the spectrum is doubled (to make the ifft real)
% --- handle PCB
Param.CTLEen = false; % enable/dis CTLE
% --- handle PCB
Normalized_Flag = input('Choose: 1 to normalized / Choose: 0 to without normalized : ');
Flag = input('Choose a PCB 1-7: '); % Choose which PCB we want to work with
switch Flag 
    case 1
        filename = 'THRU_LL_Rx1_L_Diff.mat';
        S21strc = load(filename);
        % The optimal number of Equalizer taps for SNR = 20 [dB]
         numDfeTaps = 2; 
         numFfeTaps = 16;
    case 2 
        filename = PCB_Create('10p4dB_at_26p56G_stripline.s4p','10p4dB_at_26p56G_stripline'); 
        S21strc = load(filename);
        % The optimal number of Equalizer taps for SNR = 20 [dB]
        numDfeTaps = 5; 
        numFfeTaps = 9;
    case 3
        filename = PCB_Create('15p7dB_at_26p56G_stripline.s4p','15p7dB_at_26p56G_stripline');
        S21strc = load(filename);
        % The optimal number of Equalizer taps for SNR = 20 [dB]
        numDfeTaps = 5; 
        numFfeTaps = 9;
    case 4
        filename = PCB_Create('21dB_at_26p56G_stripline.s4p','21dB_at_26p56G_stripline');
        S21strc = load(filename);
        % The optimal number of Equalizer taps for SNR = 20 [dB]
        numDfeTaps = 5; 
        numFfeTaps = 11;
    case 5
        filename = PCB_Create('26dB_at_26p56G_stripline.s4p','26dB_at_26p56G_stripline');
        S21strc = load(filename);
       % The optimal number of Equalizer taps for SNR = 20 [dB]
        numDfeTaps = 6; 
        numFfeTaps = 8;
    case 6
        filename = PCB_Create('31dB_at_26p56G_stripline.s4p','31dB_at_26p56G_stripline');
        S21strc = load(filename);
        % The optimal number of Equalizer taps for SNR = 20 [dB]
        numDfeTaps = 6; 
        numFfeTaps = 7;
    case 7 
        filename = PCB_Create('ext_thru_31dB.s4p','ext_thru_31dB');
        S21strc = load(filename);
        % The optimal number of Equalizer taps for SNR = 20 [dB]
        numDfeTaps = 1; 
        numFfeTaps = 10;
end


Param = GetFilterCoef(Param); %   Getting the filters Coefs (b and a)
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

%--- Source generation and modulation - Tx-----------------
% Generate random data source to be transmitted
[ModulatedData, Param] = dataGenerator(Param);
symbols = Param.state.symbols;
SigTransmitted = TxPS(ModulatedData,Param); % The "analog" signal 

%...............................................................
%-------- Channel and AWGN
%.............................................................
FIRlen = 50; % Fir filter length

if 1
    % SigChannel output of Tx (data)
    [SigChannel,Param] = FilterS21(SigTransmitted,S21strc,FIRlen,Param,Normalized_Flag);
else
    load(['C:\Gilad\Edi\ImpuseReponse.mat'])
    Param.ChannelResponse = h;
    FsImp =  26.56e9*8;
    channel = interp1(Param.ChannelResponse,1:FsImp/Fs:length(Param.ChannelResponse),'cubic'); channel = channel/sum(channel);
    %channel = [1 ]; % for debuging
    %freqz(channel,1,1e4,Fs)
    SigChannel = filter(channel, 1, SigTransmitted);
end

% --- Calculation of the rms value
CHANNEL_rms = rms(SigChannel);

% --- AWGN
rng(Param.Rx.NoiseSeed) % Random number generator -> rng(seed)
noise = randn(size(SigChannel))*sqrt(Fs/(2*F_cut_rx)); % *sqrt(Fs/(2*FpdCut)) - consider the Rx filter noise reduction

% --- An external function to find the optimal number of Equalizer taps for SNR = 20 [dB]
if 1 
    [Optimal_FFE_taps,Optimal_DFE_taps] = Optimize_taps(numFfeTaps_list,numDfeTaps_list,Mu,CursorPos,Param,S21strc,noise,Normalized_Flag);
end
sweepParam = SNR; % --- sweep (the SNR)

BER = nan(length(sweepParam),1);
BER_before_Dfe = nan(length(sweepParam),1);
BER_list = [];

i = 1; % FFE penalty, DFE constant
j = 1; % DFE penalty, FFE constant

for s = 1:NumOfCTLE+1
    i = 1; % FFE penalty, DFE constant
    j = 1; % DFE penalty, FFE constant
    while i <= length(numFfeTaps_list) && j <= length(numDfeTaps_list)
        for Idx = 1:length(sweepParam)
    
            % --- the channel attenuation
            SigChannel = SigChannel/rms(SigChannel); 
            SigChannel = SigChannel*10.^(+SNR(Idx)/20); 
            
            % --- receiver part
            SigRx = Rx(SigChannel,noise,Param);
    
            % --- CTLE
    %         if Idx == 0
            ImpusleResCTLE = AddCTLE_MinBER(SigRx,Param);
            CTLE = ImpusleResCTLE(:, s);
            SigRx = filter(CTLE,1,SigRx);
    %         end
    
            if 0
                eyediagram(SigRx(103:1e3),sps)
            end
    
            % --- sample at Fs = Rs
            [SigRxUsmp,~] = Sampler(SigRx,sps);
            SigRxUsmp = SigRxUsmp/max(abs(SigRxUsmp))*max(ModulatedData); % normalizing
    
            % --- FFE/DFE
            [~,w,c,Error,dataOut] = DfeLms(SigRxUsmp,ModulatedData,numFfeTaps_list(i),Mu,CursorPos,numDfeTaps_list(j));
            %     load('w_and_c')
    
            % ---using the function with a fix FFE/DFE weights
            SigDfeOut = DfeLms(SigRxUsmp,dataOut,numFfeTaps_list(i),Mu,CursorPos,numDfeTaps_list(j),0,w,c,1);
            if 0
                eyediagram(SigDfeOut(103:1e3),sps) 
            end
        
            % --- deModulation
            DataRx = deMod(SigDfeOut,Param);
            DataRx_before_Dfe = deMod(SigRxUsmp,Param);
    
            %--- BER calculation
            [BER(Idx), ~] = getBER(DataRx,Param);
            [BER_before_Dfe(Idx),~] = getBER(DataRx_before_Dfe,Param);
        
    %         figure(100);
    %         semilogy(SNR(1:Idx),BER_before_Dfe(1:Idx),SNR(1:Idx),BER(1:Idx));
    %         legend('Before DFE','After DFE')
        end
    
        BER_list = [BER_list, BER];
    
        if Checker == 1
            i = i + 1;
            if i == length(numFfeTaps_list) + 1
                j = j + 1;
                i = 1;
            end
        elseif Checker == 2
            j = j + 1;
            if j == length(numDfeTaps_list) + 1
                i = i + 1;
                j = 1;
            end
        elseif Checker == 3
            i = i + 1;
            j = j + 1;
        end
        
    end

    if (Checker == 1 || Checker == 2)
    Find_Penalty(SNR, BER_list,Checker,numFfeTaps_list,numDfeTaps_list);
    end
    
    figure(s+20)
    semilogy(SNR,BER_list(:, 1+(s-1)*6:s*6),SNR,BER_before_Dfe);
    xlim([0 30])
    title('BER performance')
    legend(['FFE taps = ',num2str(numFfeTaps_list(1)) , ',   DFE taps = ', num2str(numDfeTaps_list(1))], ...
           ['FFE taps = ',num2str(numFfeTaps_list(2)) , ',   DFE taps = ', num2str(numDfeTaps_list(2))], ...
           ['FFE taps = ',num2str(numFfeTaps_list(3)) , ',     DFE taps = ', num2str(numDfeTaps_list(3))], ...
           ['FFE taps = ',num2str(numFfeTaps_list(4)) , ',   DFE taps = ', num2str(numDfeTaps_list(4))], ...
           ['FFE taps = ',num2str(numFfeTaps_list(5)) , ',     DFE taps = ', num2str(numDfeTaps_list(5))], ...
           ['FFE taps = ',num2str(numFfeTaps_list(6)) , ',     DFE taps = ', num2str(numDfeTaps_list(6))], ...
           'BER before DFE','FontSize',13)
    ylabel('BER'); grid on
    xlabel('SNR [dB]')
    
    
    % figure(100)
    % ylabel('BER'); grid on
    % xlabel('SNR [dB]')
    
    % Plot results
    figure(100*s);
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

end

% if (Checker == 1 || Checker == 2)
%     Find_Penalty(SNR, BER_list,Checker,numFfeTaps_list,numDfeTaps_list);
% end
% 
% figure(105)
% semilogy(SNR,BER_list,SNR,BER_before_Dfe);
% xlim([0 30])
% title('BER performance')
% legend(['FFE taps = ',num2str(numFfeTaps_list(1)) , ',   DFE taps = ', num2str(numDfeTaps_list(1))], ...
%        ['FFE taps = ',num2str(numFfeTaps_list(2)) , ',   DFE taps = ', num2str(numDfeTaps_list(2))], ...
%        ['FFE taps = ',num2str(numFfeTaps_list(3)) , ',     DFE taps = ', num2str(numDfeTaps_list(3))], ...
%        ['FFE taps = ',num2str(numFfeTaps_list(4)) , ',   DFE taps = ', num2str(numDfeTaps_list(4))], ...
%        ['FFE taps = ',num2str(numFfeTaps_list(5)) , ',     DFE taps = ', num2str(numDfeTaps_list(5))], ...
%        ['FFE taps = ',num2str(numFfeTaps_list(6)) , ',     DFE taps = ', num2str(numDfeTaps_list(6))], ...
%        'BER before DFE','FontSize',13)
% ylabel('BER'); grid on
% xlabel('SNR [dB]')
% 
% 
% % figure(100)
% % ylabel('BER'); grid on
% % xlabel('SNR [dB]')
% 
% % Plot results
% figure(500);
% subplot(2,1,1)
% semilogy((Error.^2));
% title(['LMS Adaptation Learning Curve Using Mu = ', num2str(Mu)]);
% xlabel('Iteration Number');
% ylabel('Output Estimation Error in dB');
% subplot(2,1,2)
% stem(w); hold on ; stem(length(w)+1:length(w)+ length(c),c,'r') ; hold off
% title('DFE taps')
% % subplot(3,1,3)
% % stem(real(channel)); title('Channel Taps')

