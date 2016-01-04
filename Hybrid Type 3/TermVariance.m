function tv = TermVariance(Data)
% this function compute TermVariance

[S, nVar] = size(Data);

tv=zeros(nVar,1);

for  i=1:nVar
    
    tv(i) = (1/S)*sum(var(Data(:,i)));

end

tv = (tv-min(tv))./(max(tv)-min(tv));
tv(:,2) = tv(:,1);

