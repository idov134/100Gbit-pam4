function [BER, State] = getBER(symbolsDet,Param)

% symbolsTx = sim.state.dataGenerator.symbols;
% symbolsAlinged = nan(size(symbolsTx));
% symbolsDetAlinged = nan(size(symbolsTx));
symbolsTx = Param.state.symbols;
shift = nan(1,size(symbolsTx,2));
SER = nan(1,size(symbolsTx,2));
IgnoreEdge = 100;
for Idx = 1:size(symbolsTx,2)
    [symbolsDetAlinged,symbolsAlinged,shift(Idx)] = AlignDataGilad(symbolsDet(:,Idx),symbolsTx(:,Idx),1);
    SER(Idx) = sum(symbolsDetAlinged(1+IgnoreEdge:end-IgnoreEdge) ~= symbolsAlinged(1+IgnoreEdge:end-IgnoreEdge))/(size(symbolsTx,1)-2*IgnoreEdge);
end
State.SER = SER;
BER = SER/log2(Param.M); % assuming only one bit error per symbol (ueing grey code...) 
