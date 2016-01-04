function normTau = NormTau( tau )
%NORMTAU Summary of this function goes here
%   Detailed explanation goes here

    normTau =zeros(size(tau));

    for i=1:length(tau)

        normTau(i, :) = (tau(i,:))./sum(tau(i,:));

    end


end

