%% start of Program

clc;
clear;
close all;

%% Program Defination

Data  = load('Data1.mat');                      % Load Data

Data.X = NormData(Data);                       % Normalize Data

nVar=size(Data.X,2);                           % Number of Decision Variables

GlobalIt = 50;

BestCostGlobal = zeros(GlobalIt,1);

%% Problem Definition BBO

VarSize=[1 nVar];   % Decision Variables Matrix Size

VarMin=1;          % Lower Bound of Decision Variables
VarMax=2;          % Upper Bound of Decision Variables

%% BBO Parameters

BBOMaxIt=20;          % Maximum Number of Iterations

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

ACOMaxIt=10;       % Maximum Number of Iterations

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

eta=TermVariance(Data.X)+eps;      % Heuristic Information Matrix

tau=tau0*ones(nVar, 2);         % Phromone Matrix

BestCostACO=zeros(ACOMaxIt,1);       % Array to Hold Best Cost Values

%% BBO Initialization

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
    
    [pop(i).Cost,pop(i).Features]=mRMR_BBO(pop(i).Position,nVar,Data);
	
    if pop(i).Cost> BestSol.Cost
        BestSol=pop(i);
    end
end

for j = 1 : length(pop)
    
    Prob(j) = 1 / length(pop); 
    
end

% Sort Population
[~, SortOrder]=sort([pop.Cost], 'ascend');
pop=pop(SortOrder);
BBOpop = pop;

[lambda, mu] = GetLambdaMu(pop, I, E);

% Best Solution Ever Found
BestSol=pop(1);

% Array to Hold Best Costs
BestCostBBO=zeros(BBOMaxIt,1);
  tic
for iter =1:GlobalIt
    BBOpop = pop;

    %% BBO Main Loop
     for it=1:BBOMaxIt
        
        newpop = BBOpop;
        
        Prob = CalcMutRate( BBOpop, lambda, mu, Prob, dt );

            for i=1:nPop
                for k=1:length(BBOpop(i).Features)
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

            % Select Next Iteration Population
            BBOpop=[BBOpop(1:nKeep)
                 newpop(1:nNew)];

             % Get new Lambda and Mu
            [lambda, mu] = GetLambdaMu(BBOpop, I, E);

            % Sort Population
            [~, SortOrder]=sort([BBOpop.Cost], 'ascend');
            BBOpop=BBOpop(SortOrder);

            % Update Best Solution Ever Found
            BestSol=BBOpop(1);

            % Store Best Cost Ever Found
            BestCostBBO(it)=BestSol.Cost;

            % Show Iteration Information
            disp(['Iteration ' num2str(it) ': Best Cost BBO = ' num2str(BestCostBBO(it))]);
        end


    %% ACO Initialization

    % Update Phromones
    update.Costs = NormCosts(BBOpop, BestSol);

    update.Ants = [BestSol; BBOpop ];
    
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

    % Empty Ant
    empty_ant.Position=[];
    empty_ant.Features=[];
    empty_ant.Cost=[];


    % Ant Colony Matrix
    ant=repmat(empty_ant,nAnt,1);

    % Best Ant
    BestAnt.Cost=+inf;

    ACOpop=[];
    for i=1:ACOMaxIt
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

        update.Costs = NormCosts(ant, BestAnt);

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
        BestCostACO(i)=BestAnt.Cost;

        % Show Iteration Information
        disp(['Iteration ' num2str(i) ': Best Cost ACO = ' num2str(BestCostACO(i))]);
        
        ACOpop =[ACOpop
            ant];

    end

    temppop = [ACOpop; BBOpop];

    [~, SortOrder]=sort([temppop.Cost], 'ascend');
    temppop=temppop(SortOrder);


    % Select Next Iteration Population
    pop=[pop(1:nKeep)
        temppop(1:nNew)];

    

    % Sort Population
    [~, SortOrder]=sort([pop.Cost], 'ascend');
    pop=pop(SortOrder);
    
    [lambda, mu] = GetLambdaMu(pop, I, E);
    
    BestCostGlobal(iter) =  pop(1).Cost;
    
    % Show Iteration Information
    disp(['************************Iteration ' num2str(iter) ': Best Cost  = ' num2str(pop(1).Cost) '************************']);
    
    if (iter==1)
        TimeOfEachIteration = toc;
        
    end
    
end
        
    

%% Results

TimeOfEachIteration
figure;
plot (BestCostGlobal, 'LineWidth', 2);
xlabel ('Iteration');
ylabel ('Best Cost');
Out = ClassifyFunction(pop(1).Features,nVar,Data);

