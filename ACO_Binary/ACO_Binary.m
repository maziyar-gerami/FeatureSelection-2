%% Start of Program

clc
clear
close all
global Data nVar

%% Problem Definition

Data  = load('Data1.mat');                       % Load Data

Data.X = NormData(Data);                        % Normalize Data

nVar=size(Data.X,2);                            % Number of Decision Variables

CostFunction=@(x)mRMR(x,nVar,Data);             % Cost Function

%% ACO Parameters

MaxIt=20;               % Maximum Number of Iterations

nAnt=50;              % Number of Ants (Population Size)

Q=1;

q0=0.7;             % Exploitation/Exploration Decition Factor

tau0=1;             % Initial Phromone

alpha=0.7;          % Phromone Exponential Weight
beta=0.3;           % Heuristic Exponential Weight

rho=0.7;           % Evaporation Rate

%% Initialization

N=[0 1];

eta=TermVariance(Data.X)+eps;      % Heuristic Information Matrix
%eta=ones(nVar,2);

tau=tau0*ones(nVar, 2);         % Phromone Matrix

BestCost=zeros(MaxIt,1);       % Array to Hold Best Cost Values

% Empty Ant
empty_ant.Position=[];
empty_ant.Features=[];
empty_ant.Cost=[];


% Ant Colony Matrix
ant=repmat(empty_ant,nAnt,1);

% Best Ant
BestAnt.Cost=+inf;

%% ACO Main Loop
tic
for i=1:MaxIt
    %Move Ants
    for k=1:nAnt
        
        for n=1:nVar
            
           q= rand;
           
            if(q<=q0)
                
                [~, idx] = max((tau(n,:)).^alpha.*(eta(n,:)).^beta);
                
            else
                
                P = (((tau(n,:)).^alpha.*(eta(n,:)).^beta)./sum((tau(n,:)).^alpha.*(eta(n,:)).^beta));
                
                P = P/sum(P);
                
                P = P';
                
                idx = RouletteWheelSelection(P);
                
            end
            
            ant(k).Position(n) = idx;
            
        end
        
        [ant(k).Cost,ant(k).Features]=mRMR(ant(k).Position, nVar, Data);
        
    end
    
    
    [~, SortOrder]=sort([ant.Cost], 'ascend');
    ant=ant(SortOrder);
    
    if ant(1).Cost<BestAnt.Cost
        BestAnt = ant(1);
    end

    update.Costs = NormAntCosts(ant, BestAnt);
    
    update.Ants = [BestAnt; ant ];
    
    % update best path
    for j=1:nVar    
                
                tau(j, update.Ants(1).Position(j))= tau(j, update.Ants(1).Position(j))+ rho*(update.Ants(1).Cost);
                
    end
    
    % update other paths
    for k=2:20
        
        for j=1:nVar    
                
                tau(j, update.Ants(k).Position(j))= tau(j, update.Ants(k).Position(j))+ update.Ants(k).Cost;
                
        end
        
    end
     
    
    
    % Evaporation
    tau=(1-rho)*tau;
    
    % Store Best Cost
    BestCost(i)=BestAnt.Cost;
    
    % Show Iteration Information
    disp(['Iteration ' num2str(i) ': Best Cost = ' num2str(BestCost(i))]);
    
    if (i==1)
        TimeOfEachIteration = toc;
        
    end
        
end


%% Results
TimeOfEachIteration
figure;
plot (BestCost, 'LineWidth', 2);
xlabel ('Iteration');
ylabel ('Best Cost');
Out = ClassifyFunction(BestAnt.Features,nVar,Data);
