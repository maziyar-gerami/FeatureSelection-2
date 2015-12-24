function [lambda, mu] = GetLambdaMu(Population, I, E)

% Compute immigration rate and extinction rate for each species count.
% lambda(i) is the immigration rate for individual i.
% mu(i) is the extinction rate for individual i.
for i = 1 : length(Population)
    
    mu(i) = I * (1 - i / length(Population));
    
    lambda(i) = E * (i / length(Population));
    
end
end