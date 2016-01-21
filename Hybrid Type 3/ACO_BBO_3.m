%% start of Program

clc;
clear;
close all;

%% Program Defination

Data  = load('Data1.mat');                      % Load Data

Data.X = NormData(Data);                       % Normalize Data

nVar=size(Data.X,2);                           % Number of Decision Variables

GlobalIt = 5;

BGC =[];

%% Problem Definition BBO

VarSize=[1 nVar];   % Decision Variables Matrix Size

VarMin=1;          % Lower Bound of Decision Variables
VarMax=2;          % Upper Bound of Decision Variables

%% BBO Parameters

nPop=50;            % Number of Habitats (Population Size)

KeepRate=0.1;                   % Keep Rate
nKeep=round(KeepRate*nPop);     % Number of Kept Habitats

nNew=nPop-nKeep;                % Number of New Habitats

% Migration Rates
I = 1;                          % Maximum Imegration Rate
E = 1;                          % Maximum Emegration Rate
P = 12;                         % Maximum Size of Each Habitat

Mmax =0.05;

dt =1;

%% ACO Parameters

nAnt=50;        % Number of Ants (Population Size)

Q=1;

q0=0.7;         % Exploitation/Exploration Decition Factor

tau0=1;         % Initial Phromone

alpha=0.7;        % Phromone Exponential Weight
beta=0.3;         % Heuristic Exponential Weight

rho=0.7;       % Evaporation Rate



% Empty Habitat
habitat.Position=[];
habitat.Cost=[];
habitat.Features=[];

% Create Habitats Array
ACOpop=repmat(habitat,nPop,1);
BestGlobalCost =+inf;


%% BBO Initialization

% Empty Habitat
habitat.Position=[];
habitat.Cost=[];
habitat.Features=[];

% Create Habitats Array
pop=repmat(habitat,nPop,1);

% Best Habitat
BestSol.Cost=-inf;

% Initialize Habitats

%% ACO Initialization

eta=TermVariance(Data.X)+eps;      % Heuristic Information Matrix

tau=tau0*ones(nVar, 2);         % Phromone Matrix
     
tic        
for it=1:GlobalIt
        BestCostACO=zeros(GlobalIt,1);    % Array to Hold Best Cost Values

        % Empty Ant
        empty_ant.Position=[];
        empty_ant.Features=[];
        empty_ant.Cost=[];
        
        % Ant Colony Matrix
        ant=repmat(empty_ant,nAnt,1);

        % Best Ant
        BestAnt.Cost=+inf;

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

            [ant(k).Cost,ant(k).Features]=mRMR_ACO(ant(k).Position, nVar, Data);

        end


        [~, SortOrder]=sort([ant.Cost], 'ascend');
        ant=ant(SortOrder);

        if ant(1).Cost<BestAnt.Cost
            BestAnt = ant(1);
        end
        
        ACOpop =ant;
        

        % Store Best Cost
        BestCostACO(it)=BestAnt.Cost;

        % Show Iteration Information
        disp(['Iteration ' num2str(it) ': Best Cost ACO= ' num2str(BestCostACO(it))]);       


        for j = 1 : length(ACOpop)

            Prob(j) = 1 / length(ACOpop); 

        end

        % Sort Population
        [~, SortOrder]=sort([ACOpop.Cost], 'ascend');
        ACOpop=ACOpop(SortOrder);
        ACOpop = ACOpop (1:50,:);

        [lambda, mu] = GetLambdaMu(pop, I, E);

        % Array to Hold Best Costs
        BestCostBBO=zeros(GlobalIt,1);

        %% BBO Main Loop

        
        newpop=ACOpop;
        
        Prob = CalcMutRate( newpop, lambda, mu, Prob, dt );
        
        for i=1:nPop
                for k=1:length(ACOpop(i).Features)
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
                if (sum( newpop(i).Position))==0
                    newpop(i).Cost =+inf;

                % Evaluation
                else
                    [newpop(i).Cost,newpop(i).Features]=mRMR_BBO(newpop(i).Position,nVar,Data);
                end

        end
        
       % Sort New Population
       [~, SortOrder]=sort([newpop.Cost],'ascend');
        newpop=newpop(SortOrder); 
        
        BestCostBBO = newpop(1).Cost;
        
        disp(['Iteration ' num2str(it) ': Best Cost BBO= ' num2str(BestCostBBO)]); 
        
        BBOpop = newpop;
        
        % Select Next Iteration Population
        pop=[ACOpop(1:nKeep)
            BBOpop(1:nNew)];
        

        % Sort Population
        [~, SortOrder]=sort([pop.Cost], 'ascend');
        pop=pop(SortOrder);
        
        BestSol = pop(1);

      
        BestCostGlobal(it) =  BestSol.Cost;

        update.Costs = NormCosts(pop, BestSol);

        update.Ants = [BestSol; pop];
        
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

        
        if  BestSol.Cost < BestGlobalCost
            BestGlobalCost = BestSol.Cost; 
            BestGlobalSol = BestSol;
        end
        
        BGC = [BGC, BestGlobalCost];

        disp(['************* Iteration ' num2str(it) ': Best Cost  = ' num2str(BestGlobalCost)]);
        
        if (it==1)
            
            TimeOfEachIteration = toc;
        
        end
end

Time = toc;

%% Results
TimeOfEachIteration
figure;
plot (BGC, 'LineWidth', 2);
xlabel ('Iteration');
ylabel ('Best Cost');
Out = ClassifyFunction(pop(1).Features,nVar,Data);

