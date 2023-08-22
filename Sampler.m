function [Out,IndexMax] = Sampler(In,sps)

temp = reshape(In(1:length(In) - rem(length(In),sps)),sps,[]);
[~,IndexMax] = max(std(temp'));
Out = In(IndexMax:sps:end);
end

