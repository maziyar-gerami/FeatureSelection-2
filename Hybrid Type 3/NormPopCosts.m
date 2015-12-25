function Costs = NormPopCosts( Pop, BestAnt )


    
    
    Costs=[];

    for i=1:length(Pop)
       
        Costs = [Costs; Pop(i).Cost];
        
    end
    
    
    normData = max(Costs) - min(Costs);               % this is a vector
    
    normData = repmat(normData, [length(Costs) 1]);    % this makes it a matrix
    
    minmatrix= repmat(min(Costs), [length(Costs) 1]); % of the same size as A
                                                         
    Costs = (Costs-minmatrix)./normData;         % your normalized matrix

    
    
    
end

