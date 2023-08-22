function [Out,w,c,Error,Data,IndexCnvrgnc] = DfeLms(Sig,Data,numFfeTaps,Mu,CursorPos,numDfeTaps,plotEn,w,c,UpdateEn)

[Sig,Data,idxOut] = AlignData(Sig,Data,1);
Sig = Sig(CursorPos:end);

numPoints = min([length(Data), length(Sig)]);
N = numPoints - numFfeTaps;
Out = nan(N,1);
Error = nan(N,1);

if nargin < 7
    plotEn = 0;
end

if nargin < 8 
    w = zeros(numFfeTaps+1,1);% + i*zeros(numTaps+1,1); % FFE weights
    c = zeros(numDfeTaps,1); % DFE weights
    UpdateEn = 1;
end


M = max([numFfeTaps,numDfeTaps]);
% LMS Adaptation
for n  = M+1 : numPoints
    
    % select part of training input
    in = Sig(n : -1 : n-numFfeTaps) ;
    DataTmp = Data(n-1:-1:n-numDfeTaps);
    
    Out(n) = w'*in + c'*DataTmp;
    
    % compute error
    Error(n) = Data(n)-Out(n);
    
    % update taps
    if UpdateEn 
        w = w + Mu*(real(Error(n)*conj(in)));% - 1i*imag(Error(n)*conj(in)) );
        c = c + Mu*(Error(n)*DataTmp);%
    end
    if plotEn
        figure(555)
        subplot(3,1,1)
        plot(Error(M+1:n))
        subplot(3,1,2)
        stem(w)
        subplot(3,1,3)
        stem(c)
    end
    
end
Out = Out(M+1:numPoints);
Data = Data(M+1:numPoints);
% Find convergence point
ErrorMA = filter(ones(floor(length(Error)/100),1),1,Error.^2);
IndexCnvrgnc = find(ErrorMA < min(ErrorMA)*1.1,1,'first');

end

