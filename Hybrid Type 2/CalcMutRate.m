function Prob = CalcMutRate( Population, lambda, mu, Prob, dt )

        ProbDot(1) = -(lambda(1)+mu(1))*Prob(1)+ mu(2)*Prob(2);

        for j = 2 : length(Population)-1
            % Compute lambda for one less than the species count of habitat i.
            lambdaMinus = lambda(j-1);
            % Compute mu for one more than the species count of habitat i.
            muPlus = mu(j+1);
            % Compute Prob for one less than and one more than the species count of habitat i.
            % Note that species counts are arranged in an order opposite to that presented in
            % MacArthur and Wilson's book - that is, the most fit
            % habitat has index 1, which has the highest species count.
            if j < length(Population)
                ProbMinus = Prob(j+1);
            else
                ProbMinus = 0;
            end
            if j > 1
                ProbPlus = Prob(j-1);
            else
                ProbPlus = 0;
            end
            ProbDot(j) = -(lambda(j) + mu(j)) * Prob(j) + lambdaMinus * ProbMinus + muPlus * ProbPlus;
        end
        
        ProbDot(length(Population)) = -(lambda(length(Population)-1)+mu(length(Population)-1))*Prob(length(Population)-1)+ mu(length(Population)-1)*Prob(length(Population)-1);
        % Compute the new probabilities for each species count.
        Prob = Prob + ProbDot ;
        Prob = max(Prob, 0);
        Prob = Prob / sum(Prob); 
        

end


