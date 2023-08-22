clear; close all;
PCB = sparameters('31dB_at_26p56G_stripline.s4p'); % we use DAC s21 as the difference to ADC is small and need to be changed anyway
temp = s2sdd(PCB.Parameters,1);
SparamDiff = sparameters(temp,PCB.Frequencies, 2*PCB.Impedance);
% rfplot(SparamDiff)
S21 = squeeze(SparamDiff.Parameters(1,2,:));
F = SparamDiff.Frequencies
plot(F,db(S21)); grid
save('31dB_at_26p56G_stripline','SparamDiff')
%%
FirLen = 30;
[FIR] = S21toFIR(S21,F,FirLen)