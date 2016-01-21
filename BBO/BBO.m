tic
clc;
clear;
close all;

%% Problem Definition

Data  = load('Data.mat');                      % Load Data

Data.X = NormData(Data);                       % Normalize Data

nVar=size(Data.X,2);                           % Number of Decision Variables

VarSize=[1 nVar];   % Decision Variables Matrix Size

VarMin=1;          % Lower Bound of Decision Variables
VarMax=2;          % Upper Bound of Decision Variables

%% BBO Parameters

MaxIt=50;          % Maximum Number of Iterations

nPop=50;            % Number of Habitats (Population Size)

KeepRate=0.2;                   % Keep Rate
nKeep=round(KeepRate*nPop);     % Number of Kept Habitats

nNew=nPop-nKeep;                % Number of New Habitats

% Migration Rates


I = 1;                          % Maximum Imegration Rate
E = 1;                          % Maximum Emegration Rate
P = 12;                         % Maximum Size of Each Habitat

alpha=1;

pMutation=0.1;

Mmax=1;

dt =1;
%% Initialization

% Empty Habitat
habitat.Position=[];
habitat.Cost=[];
habitat.Features=[];

% Create Habitats Array
pop=repmat(habitat,nPop,1);

% Best Habitat
BestSol.Cost=+inf;

Prob =zeros(1,nPop);

% Initialize Habitats
for i=1:nPop
    pop(i).Position=randi([VarMin VarMax],VarSize);
    
    [pop(i).Cost,pop(i).Features]=mRMR(pop(i).Position,nVar,Data);
	
    if pop(i).Cost> BestSol.Cost
        BestSol=pop(i);
    end
end

% Sort Population
[~, SortOrder]=sort([pop.Cost], 'ascend');
pop=pop(SortOrder);

% Get Lambda and Mu
[lambda, mu] = GetLambdaMu(pop, I, E);

% Best Solution Ever Found
BestSol=pop(1);

% Array to Hold Best Costs
BestCost=zeros(MaxIt,1);

% Calculating Probaiity for Intial Solutions
for j = 1 : nPop
    
    Prob(j) = 1 / length(pop); 
    
end

%% BBO Main Loop

for it=1:MaxIt
    
    newpop=pop;
    Prob = CalcMutRate( pop, lambda, mu, Prob, dt );
        
        for i=1:nPop
            for k=1:length(pop(i).Features)
                % Migration
                if rand<=lambda(i)
                    % Emmigration Probabilities
                    EP=mu;
                    EP(i)=0;
                    EP=EP/sum(EP);

                    % Select Source Habitat                  
                    flag = true;
                    while(flag)
                        
                        j=RouletteWheelSelection(EP);
                        
                        if ~(isempty(newpop(j).Features))
                           flag =false; 
                        end
                    end
                    
                    
                    
                    idx = randi((numel(newpop(j).Features)));

                    % Migration
                    newpop(i).Position(newpop(i).Features(k)) = 1;
                    
                    newpop(i).Position(newpop(j).Features(idx)) = 2;
                    
                    EmigratedFeature = newpop(j).Features(idx);
                    
                    newpop(i).Features(k) = EmigratedFeature;


                end

                % Mutation                
                Pmax = max(Prob);
                MutationRate = Mmax * (1 - Prob / Pmax);
                
                if rand<MutationRate(i)
                        if (newpop(i).Position(k))==1
                            newpop(i).Position(k)=2;
                        else
                            newpop(i).Position(k)=1;
                        end
                                    
                    end
            end

            % Is Feasible?
             
            % Evaluation

            [newpop(i).Cost,newpop(i).Features]=mRMR(newpop(i).Position,nVar,Data);

        end

        % Sort New Population
        [~, SortOrder]=sort([newpop.Cost],'ascend');
        newpop=newpop(SortOrder);

        % Select Next Iteration Population
        pop=[pop(1:nKeep)
             newpop(1:nNew)];
        

        % Sort Population
        [~, SortOrder]=sort([pop.Cost], 'ascend');
        pop=pop(SortOrder);
        
        % Get new Lambda and Mu
        [lambda, mu] = GetLambdaMu(pop, I, E);

        % Update Best Solution Ever Found
        BestSol=pop(1);

        % Store Best Cost Ever Found
        BestCost(it)=BestSol.Cost;

        % Show Iteration Information
        disp(['Iteration ' num2str(it) ': Best Cost BBO = ' num2str(BestCost(it))]);
    
end
Time = toc;

%% Results
Time
figure;
plot (BestCost, 'LineWidth', 2);
xlabel ('Iteration');
ylabel ('Best Cost');
Out = ClassifyFunction(BestSol.Features,nVar,Data);




