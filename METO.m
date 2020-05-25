%% METO Algorithm

function[Sol, SolutionConverganceCurve, Feval] = METO(CostFunction,nVar, VarMin, VarMax, vType, nPop, iteration, max_f_evaluation)
%  tic   
 %% lower and upper bound

%% Selecting the number of  bits to encode a variable
if length(VarMin) < 2 % Upper and Lower Bound on nVar variables
   VarMin = repmat(VarMin,1,nVar); 
   VarMax = repmat(VarMax,1,nVar); 
end

TypeReal = find(vType ==1);
TypeInteger =find(vType ==2);
TypeBinary = find(vType ==3);

bits = zeros(1,nVar); % number of bits to represent each variable
for i = 1:nVar
   if sum(i == TypeReal)>0
       bits(i) =ceil(log2((VarMax(i) - VarMin(i))*10000));
   elseif sum(i == TypeInteger)>0
       bits(i) =ceil(log2(VarMax(i) - VarMin(i)));
   elseif sum(i == TypeBinary)>0
      bits(i) =1;  
   end  
end

% for i = 1 : length(TypeReal)  
%   bits(i) =ceil(log2((VarMax(i) - VarMin(i))*10000)); 
% end 
% for i = 1 : length(TypeInteger)  
%   bits(i) =ceil(log2(VarMax(i) - VarMin(i))); 
% end 
% for i = 1 : length(TypeBinary)  
%   bits(i) =1; 
% end 

%% Defining Species 
species =2; % number of species should be even number such as 2,4,6,8,...
% Generally for 2 species METO gives good result. It can be increased based
% on the problem complexity.
Offsprings = 2; % Number of F1 Generation offspring from a parent pair
nPop = species* round(nPop/species); % making population size feasible to distribute the candidated in the species.
p =nPop/species;% Number of Strands  
[temp_R, temp_V, temp_P, temp_f, temp_X] = traits(); % Initialize variables
f_eval =0; % initialize function evaluation number

[pop] = InitialPopulation(species, p,nVar, ...
      bits, temp_R, temp_P, temp_f, temp_X, temp_V, CostFunction, VarMin, VarMax); % Initialization of species population
  
SolutionConverganceCurve = min(pop.(temp_P{1}).L_best.f); % Initialization of Best Solution

Feval  = [];
for iter = 1: iteration % How many Evolution Epochs you want to run the METO algorithm 
[polination_ID, ~] = polination (species); % Select the species for Cross-Breeding
mut_G_a = 1./bits; % Lower limit of Epimutation probability
% Heredity will transfer if it is better than the parents genes
a_Mp = 0.94; % Lower limit of Mendelian Probability for heredity transfer
b_Mp = 0.99; % Upper limit of Mendelian Probability for heredity transfer
for i=1:species
% Parameters of METO
F1_offspring_rate = max(0.3,rand);% Cross Breeding rate (CBR) (It decides that how many parents are available for CBR)
F2_offspring_rate = max(0.5,rand);% Self Breeding rate (SBR) (It decides that how many parents are available for SBR)
if nVar >= 10 % This condition will tru if number of variables are greater than 9
tau = 1/nVar + (randperm(round(nVar/4),1)/nVar)*rand; % variables selection rate for mutation
selectVariable = 1/nVar + (randperm(round(nVar/4),1)/nVar)*rand; %rand; % Variable Selection Rate for AS strand formation
else
tau = max(0.1,rand); % for less number of variables
selectVariable = rand; % for less number of variables
end
rr = rand; % Number of bits to flip in a segment of SS to form ASS using FLipper Operator
mut_G_b =  max(1./bits, 0.01*rand + (1 - 0.91*(f_eval/max_f_evaluation)^0.005)); % Upper limit of Epimutation probability

[~,IdBest] =sort(pop.(temp_P{polination_ID(i,1)}).R.(temp_f{polination_ID(i,1)}));
ChfunNumF1 = sort(IdBest(1:round(p*F1_offspring_rate))); % Selected ID of Plants for Cross Breeding in the species populaiton
% In the both breeder species, Plants fro Cross Breeding are at the same position.
%% Taking best recessive genes from the two strands of DNA

id_L = pop.(temp_P{polination_ID(i,2)}).V.(temp_f{polination_ID(i,2)})(ChfunNumF1) <= pop.(temp_P{polination_ID(i,1)}).L_best.f(ChfunNumF1);
id_L = ChfunNumF1(id_L); % Id for updating the heredity if genes of ASS of breeder is better 

if sum(id_L) > 0  % As a Cross Breeding between two species parents strands good gens are presereved by nature     
%Phenotype/ Binary inheretence in the DNA strand for ith species
pop.(temp_P{polination_ID(i,1)}).L_best.string(id_L,:) = pop.(temp_P{polination_ID(i,2)}).V.(temp_V{polination_ID(i,2)})(id_L,:); 
% Evaluted Inherent DNA strand for ith species
pop.(temp_P{polination_ID(i,1)}).L_best.f(id_L,:)  =  pop.(temp_P{polination_ID(i,2)}).V.(temp_f{polination_ID(i,2)})(id_L,:);                
% Genotype/ Real valued inheretence in the DNA strand for ith species
pop.(temp_P{polination_ID(i,1)}).L_best.X (id_L,:)= pop.(temp_P{polination_ID(i,2)}).V.(temp_X{polination_ID(i,2)})(id_L,:);  
end

%% Mendelian Process: F1 and F2 generation offspring and heredity transfer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% F1-generation Offsprings ("F1MendelianOffspring" function produces multiple F1 offspring)
[offspring_F1] = F1MendelianOffspring...
                            (pop.(temp_P{polination_ID(i,1)}).R.(temp_R{polination_ID(i,1)})(ChfunNumF1,:),...
                             pop.(temp_P{polination_ID(i,2)}).V.(temp_V{polination_ID(i,2)})(ChfunNumF1,:),...
                             nVar, bits, length(ChfunNumF1), Offsprings);                      

for os = 1:Offsprings % Evalute all set of F1 generation offspring
[ F1_f(os).offspring,F1_X(os).offspring ]= variable_value(length(ChfunNumF1),...
                offspring_F1(os).offspring, bits, nVar, CostFunction, VarMin, VarMax);
%"variable_value" is a function to evalute genotype strands in to the
%phenotype 
 f_eval = f_eval + length(ChfunNumF1); % Number of function evaluation until this point
end

if Offsprings > 1 % taking best offspring from all sets of offspring
    for F1offspring =2:Offsprings 
        bestOffspring =  F1_f(1).offspring >  F1_f(F1offspring).offspring; % Binary coded offspring
        F1_X(1).offspring(bestOffspring,:) = F1_X(F1offspring).offspring(bestOffspring,:);% Real offspring
        F1_f(1).offspring(bestOffspring,:) =F1_f(F1offspring).offspring(bestOffspring,:); % Evaluation of offspring
        offspring_F1(1).offspring(bestOffspring,:) = offspring_F1(F1offspring).offspring(bestOffspring,:);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% replace bad parents by good F1 generation offspring
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Id = F1_f(1).offspring <= pop.(temp_P{i}).R.(temp_f{i})(ChfunNumF1);% Select Id of better offspring
IdchangeF1 = ChfunNumF1(Id);% Select Better F1 offspring 
if sum(IdchangeF1 ) > 0
pop.(temp_P{i}).R.(temp_R{i})(IdchangeF1,:) = offspring_F1(1).offspring(Id,:);
pop.(temp_P{i}).R.(temp_f{i})(IdchangeF1,:) = F1_f(1).offspring(Id,:);
pop.(temp_P{i}).R.(temp_X{i})(IdchangeF1,:) = F1_X(1).offspring(Id,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make AS-strand (ASS) for F1 generation offspring
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,Idworst] = sort(F1_f(1).offspring); % This mechanism will try to evolve bad plants through ASS
worst = round(length(Idworst)*max(0.2,rand)); % Few worst plants are selected for evolution, 
%They will Breed with next Evolution population of SS
Idworst = sort(Idworst(1:worst));
IdchangeF1 = ChfunNumF1(Idworst);

if ~isempty(IdchangeF1) > 0 % "ASstrand" function convert SS (3'-5') to ASS (5'-3')

 V_MC_F1 = offspring_F1(1).offspring(Idworst,:);
[V_MC_F1, F1tempVf, F1tempVX, feval] = ASstrand(length(IdchangeF1), V_MC_F1, offspring_F1(1).offspring(Idworst,:),...
                                         nVar, bits,CostFunction, VarMin, VarMax,  rr, selectVariable);
f_eval=f_eval +feval; % Number of function evaluation until this point

id = F1tempVf <= pop.(temp_P{i}).L_best.f(IdchangeF1); % replace heredity of offspring through ASS
IdchangeF1 = ChfunNumF1(id);
if sum(id) > 0    
pop.(temp_P{i}).L_best.string(IdchangeF1,:) = V_MC_F1(id,:); 
pop.(temp_P{i}).L_best.f(IdchangeF1,:) =  F1tempVf(id,:);                
pop.(temp_P{i}).L_best.X (IdchangeF1,:)= F1tempVX(id,:);
end
end
 
id_L = pop.(temp_P{i}).R.(temp_f{i}) <= pop.(temp_P{i}).L_best.f;% change heredity through F1-generation offspring
if sum(id_L) > 0   
pop.(temp_P{i}).L_best.string(id_L,:) = pop.(temp_P{i}).R.(temp_R{i})(id_L,:); 
pop.(temp_P{i}).L_best.f(id_L,:) =  pop.(temp_P{i}).R.(temp_f{i}) (id_L,:);                
pop.(temp_P{i}).L_best.X (id_L,:)= pop.(temp_P{i}).R.(temp_X{i})  (id_L,:);
end

%% Select good F1 offspring to generate F2 offspring

[xBreeded] =pop.(temp_P{i}).R.(temp_f{i})(ChfunNumF1);
[~,IdBest] = sort(xBreeded);
select = sort(IdBest(1:round(length(ChfunNumF1)*F2_offspring_rate)));
ChfunNumF2 = ChfunNumF1(select); % Id from F1 generation offspring which will produce F2 generation offspring
Epirate =ChfunNumF2; % Only those plants are subjected to Epimutate which are going to produce F2-offspring 
% "Epirate" is the Epimutation rate

%% Epimutation in Heredity before geneating F2 generation offspring:

 if rand > 0 % This can be the topic of research
[pop.(temp_P{i}).L_best,feval] = Epimutation(nVar, bits, mut_G_a, mut_G_b , Epirate,...
                                pop.(temp_P{i}).L_best, CostFunction,VarMin, VarMax, tau);
f_eval = f_eval + feval;  % Number of function evaluation until this point 
 end
 
%% F2-generation Offsprings ("F2MendelianOffspring" function produces F2 offspring)
% Multiple F2 offspring generation can improve the optimizer and the topic
% of research

     selectedF2fromF1 = ChfunNumF2; % It will be used when replacing population 

     SelfBreedOrganism = pop.(temp_P{i}).R.(temp_R{i})(selectedF2fromF1,:); % Selected F1 offspring
     SelfBreedOrganismFitness =   pop.(temp_P{i}).R.(temp_f{i})(selectedF2fromF1,:)  ; 
     Heredity = pop.(temp_P{i}).L_best.string(selectedF2fromF1,:); % Selected heredity corresponding to selected F1 offspring
     HeredityFitness = pop.(temp_P{i}).L_best.f(selectedF2fromF1,:);
     
[MC_F2, Index, feval] = F2MendelianOffspring... % This function produces F2 generation offspring
         (SelfBreedOrganism, Heredity,SelfBreedOrganismFitness,HeredityFitness, a_Mp,b_Mp);
   
if sum( Index)<1 % This is for avoiding error with null offspring
  MC_F2_f =[];
 MC_F2_X =[]; 
else
[MC_F2_f, MC_F2_X] ...
         = variable_value(size(MC_F2,1), MC_F2, bits, nVar, CostFunction, VarMin, VarMax);
f_eval = f_eval + feval;  % Number of function evaluation until this point 
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% replace bad parents by good F2 generation offspring
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Id  = selectedF2fromF1(Index);
IdchangeF2 = MC_F2_f <= pop.(temp_P{i}).R.(temp_f{i})(Id); % Replacing bad parents with good F2 offspring
if sum(IdchangeF2 ) > 0
pop.(temp_P{i}).R.(temp_R{i})(Id(IdchangeF2),:) = MC_F2(IdchangeF2,:);
pop.(temp_P{i}).R.(temp_f{i})(Id(IdchangeF2),:) = MC_F2_f(IdchangeF2,:);
pop.(temp_P{i}).R.(temp_X{i})(Id(IdchangeF2),:) = MC_F2_X(IdchangeF2,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make AS-strand (ASS) from  F2 generation offspring
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,Idworst] = sort(MC_F2_f);
if length(Idworst) < 2
 worst =   1:length(Idworst); 
else
worst = max(1, round(length(Idworst)*max(0.5,rand))); % Atleast one plant is selected
end
Idworst = sort(Idworst(worst));
Idchange = Id(Idworst);

V_MC_F2 = MC_F2(Idworst,:);
 if length( Idchange) > 0
[V_MC_F2, F2tempVf, F2tempVX, feval] = ASstrand(length(Idchange), V_MC_F2, MC_F2(Idworst,:) ,...
                                         nVar, bits,CostFunction, VarMin, VarMax, rr, selectVariable);   
f_eval=f_eval +feval;  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change heredity by AS strands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
Idch = F2tempVf <= pop.(temp_P{i}).L_best.f(Idchange,:);
if sum(Idch) > 0    
pop.(temp_P{i}).L_best.string((Idchange(Idch)),:) = V_MC_F2(Idch,:); 
pop.(temp_P{i}).L_best.f((Idchange(Idch)),:) =  F2tempVf(Idch,:);                
pop.(temp_P{i}).L_best.X ((Idchange(Idch)),:)= F2tempVX(Idch,:);
end
else
     V_MC_F2 =[];
     F2tempVf =[]; 
     F2tempVX =[];
end 
IdchangeF2 = MC_F2_f <= pop.(temp_P{i}).L_best.f(Id);
if sum(IdchangeF2 ) > 0   % Replacing heredity genes by good genes of F2 generation
pop.(temp_P{i}).L_best.string(Id(IdchangeF2),:) = MC_F2(IdchangeF2,:); 
pop.(temp_P{i}).L_best.f(Id(IdchangeF2),:) =  MC_F2_f(IdchangeF2,:);                
pop.(temp_P{i}).L_best.X (Id(IdchangeF2),:)= MC_F2_X(IdchangeF2,:);
end                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Select good ASS strands from all available to form virtual ASS population
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempVf = [ pop.(temp_P{i}).V.(temp_f{i}) ; F1tempVf; F2tempVf];
V = [pop.(temp_P{i}).V.(temp_V{i}) ; V_MC_F1;V_MC_F2];
tempVX = [ pop.(temp_P{i}).V.(temp_X{i}) ; F1tempVX; F2tempVX];
[~, ida] = sort(tempVf);

 ida = sort(ida(1:p));
     pop.(temp_P{i}).V.(temp_V{i})=V(ida,:);
     pop.(temp_P{i}).V.(temp_f{i})=tempVf(ida,:);
     pop.(temp_P{i}).V.(temp_X{i})=tempVX(ida,:);


%% Shuffling the candidates in the species for removing any biasing effect [Preserving the heredity]
if rand > 0 %This condition is the research topic 
   swap = randperm(p,p);% SS, ASS and Heredity of a particular plant should be at the same position
    pop.(temp_P{i}).L_best.string = pop.(temp_P{i}).L_best.string(swap,:); 
    pop.(temp_P{i}).L_best.f = pop.(temp_P{i}).L_best.f(swap,:); 
    pop.(temp_P{i}).L_best.X = pop.(temp_P{i}).L_best.X(swap,:); 
    
    pop.(temp_P{i}).R.(temp_R{i}) = pop.(temp_P{i}).R.(temp_R{i})(swap,:); 
    pop.(temp_P{i}).R.(temp_f{i}) = pop.(temp_P{i}).R.(temp_f{i})(swap,:); 
    pop.(temp_P{i}).R.(temp_X{i}) =  pop.(temp_P{i}).R.(temp_X{i})(swap,:); 

   pop.(temp_P{i}).V.(temp_V{i}) =  pop.(temp_P{i}).V.(temp_V{i})(swap,:); 
   pop.(temp_P{i}).V.(temp_f{i}) = pop.(temp_P{i}).V.(temp_f{i})(swap,:); 
   pop.(temp_P{i}).V.(temp_X{i}) =  pop.(temp_P{i}).V.(temp_X{i})(swap,:); 
end

%%
  [a_temp(i), id_a(i)]  =  min(pop.(temp_P{i}).L_best.f ); % preserving plant with best heredity from each species
 
end

 [F_test, sol_Id] = min(a_temp); %Preserving the Pseudo-best solution from all species
if iter == 1
        Sol.str = pop.(temp_P{sol_Id}).L_best.string(id_a(sol_Id),:);
        Sol.f = pop.(temp_P{sol_Id}).L_best.f(id_a(sol_Id),:);
        Sol.X = pop.(temp_P{sol_Id}).L_best.X(id_a(sol_Id),:);
        SolutionConverganceCurve(iter+1) = Sol.f;
else
 if F_test <= Sol.f      
         Sol.str = pop.(temp_P{sol_Id}).L_best.string(id_a(sol_Id),:);
        Sol.f  = pop.(temp_P{sol_Id}).L_best.f(id_a(sol_Id),:);
        Sol.X= pop.(temp_P{sol_Id}).L_best.X(id_a(sol_Id),:);
        SolutionConverganceCurve(iter+1) = Sol.f;
else
%       Sol.f(iter)  = Sol.f(iter-1);
      SolutionConverganceCurve(iter+1) = SolutionConverganceCurve(iter);
 end
end
convergence(iter,SolutionConverganceCurve,f_eval) 
  if f_eval >= max_f_evaluation 
         break
  end 
end
% toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = convergence(iter,f,f_eval)  
% If you will comment it then no display of the iterations
disp(['Epochs =' num2str(iter) ': F_eval =' num2str(f_eval) ': Best Cost = ' num2str(f(iter))]);
end
%%
function [polination_ID, Polination] = polination (species)
% This function is for selecting the species for Cross Breeding
r=1;
temp = 1 : species;
polination_ID = nan(species,2); 
Polination = nan(species/2,2);
while r <= species/2
   random_num=round(2+(length(temp)-2)*rand);
   Polination(r,:)= [temp(1) temp(random_num)];
   temp = setdiff(temp , Polination(r,:));
   r=r+1;
end
polination_ID(Polination(:,1),:)= Polination;
polination_ID(Polination(:,2),:)= fliplr(Polination);
end

 function [temp_R, temp_V, temp_P, temp_f, temp_X] = traits()
 % Maximum 30 species can generate
temp_R = {'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9', 'R10', 'R11', 'R12', 'R13', 'R14', 'R15', 'R16',...
        'R17', 'R18', 'R19', 'R20', 'R21', 'R22', 'R23', 'R24', 'R25', 'R26', 'R27', 'R28', 'R29', 'R30'};
temp_V = {'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'V10', 'V11', 'V12', 'V13', 'V14', 'V15', 'V16',...
        'V17', 'V18', 'V19', 'V20', 'V21', 'V22', 'V23', 'V24', 'V25', 'V26', 'V27', 'V28', 'V29', 'V30'};
temp_P = {'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15', 'P16',...
        'P17', 'P18', 'P19', 'P20', 'P21', 'P22', 'P23', 'P24', 'P25', 'P26', 'P27', 'P28', 'P29', 'P30'};
temp_f = {'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'f11', 'f12', 'f13', 'f14', 'f15', 'f16', ...
        'f17', 'f18', 'f19', 'f20', 'f21', 'f22', 'f23', 'f24', 'f25', 'f26', 'f27', 'f28', 'f29', 'f30'};
temp_X = {'X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7', 'X8', 'X9', 'X10', 'X11', 'X12', 'X13', 'X14', 'X15', 'X16',...
        'X17', 'X18', 'X19', 'X20', 'X21', 'X22', 'X23', 'X24', 'X25', 'X26', 'X27', 'X28', 'X29', 'X30'};  
 end

 
 function[f, X]= variable_value(nPop, pop, bits, nVar, CostFunction, VarMin, VarMax)
 f = zeros(nPop,1); % Initialize vector having function evaluations for each plant in the species 
Var_value_pop=nan(nPop,nVar); % prelocation of the variables
X =  Var_value_pop;     
        for i = 1 : nVar  
            if bits(i)  == 1   
             X(:,i) = pop(:,i);  
            else
        Power = repmat(2.^fliplr(0:(bits(i)-1)),nPop,1); % 2^(n-1)
        Var_value_pop(:,i) = sum(pop(:,((i - 1) * bits(i)) + 1 : i*bits(i)).* Power, 2); 
        X(:,i) = repmat(VarMin(i),nPop,1)+ repmat((VarMax(i)-VarMin(i))./(2.^bits(i)-1),nPop,1).*Var_value_pop(:,i);
            end
        end 
%   X = repmat(VarMin,nPop,1)+ repmat((VarMax-VarMin)./(2.^bits-1),nPop,1).*Var_value_pop;

for i =1: size(X,1)
 f(i,:)=CostFunction(X(i,:));
end
 end


 function [pop] = InitialPopulation(species, nPop, nVar, bits,...
    temp_R, temp_P, temp_f, temp_X,temp_V, CostFunction, VarMin, VarMax)
% Initialization of SS and ASS populations for each species of plants 

 for i = 1:species    
     pop.(temp_P{i}).R.(temp_R{i}) =[];
    for v = 1:nVar
    pop.(temp_P{i}).R.(temp_R{i}) = [pop.(temp_P{i}).R.(temp_R{i}) round(rand(nPop,bits(i)))]; % random generation of the SS population 
    end
   [pop.(temp_P{i}).R.(temp_f{i}), pop.(temp_P{i}).R.(temp_X{i})] ...
       = variable_value(nPop, pop.(temp_P{i}).R.(temp_R{i}), bits, nVar, CostFunction, VarMin, VarMax); % Fitness and point calculation

pop.(temp_P{i}).L_best.string  = pop.(temp_P{i}).R.(temp_R{i}); %Initialize Heredity Genes 
pop.(temp_P{i}).L_best.f  =  pop.(temp_P{i}).R.(temp_f{i});   %Evalute Heredity Genes             
pop.(temp_P{i}).L_best.X = pop.(temp_P{i}).R.(temp_X{i}); %Convert Heredity Genome to Phenotype
                          
pop.(temp_P{i}).V.(temp_V{i})= fliplr(pop.(temp_P{i}).R.(temp_R{i}));% ASS of the initial SS population
[pop.(temp_P{i}).V.(temp_f{i}),  pop.(temp_P{i}).V.(temp_X{i})]...
       = variable_value(nPop, pop.(temp_P{i}).V.(temp_V{i}), bits, nVar, CostFunction, VarMin, VarMax);% Fitness calculation
% id_L = pop.(temp_P{i}).V.(temp_f{i}) < pop.(temp_P{i}).L_best.f;
% if ~isempty(id_L)  
% pop.(temp_P{i}).L_best.string(id_L,:)= pop.(temp_P{i}).V.(temp_V{i})(id_L,:); 
% pop.(temp_P{i}).L_best.f(id_L,:) =  pop.(temp_P{i}).V.(temp_f{i}) (id_L,:);                
% pop.(temp_P{i}).L_best.X (id_L,:) = pop.(temp_P{i}).V.(temp_X{i})  (id_L,:);   
% end
 end
 end
 
 %% F1 Generation MendelianOffspring:
function [MF1] = F1MendelianOffspring(principal_popu,  popu_replica, nVar, bits, nPop, Offsprings)
r = false(nPop,sum(bits)); 
for j = 1:nPop
   SelectDna = sort(randperm(nVar, ceil(nVar*rand)));
   for i = SelectDna
       if i ==1
       Bits = 1:bits(1);    
       else
       Bits = sum(bits(1:i-1))+1:sum(bits(1:i));
       end

       r(j,Bits) = true;
   end
end
r_not = ~r; 
for k = 1 : Offsprings
    children_F1 = principal_popu;
    if k == 1 % For generating only one offspring
        children_F1(r) = popu_replica (r);
    elseif k ==2 % For generating second offspring
        children_F1(r_not) = popu_replica (r_not);
    else
        r = false(nPop,sum(bits));
        for j = 1:nPop
            SelectDna = randperm(nVar, ceil(nVar*rand));
             for i = SelectDna
                 if i ==1
                    Bits = 1:bits(1);    
                 else
                    Bits = sum(bits(1:i-1))+1:sum(bits(1:i));
                end
                r(j,Bits) = true;
             end
        end    
    end
    MF1(k).offspring = children_F1;
end
end

%% F2 Generation MendelianOffspring:
function [F1_offspring, Index1, feval] = F2MendelianOffspring...
    (F1_offspring,  heredity, F1_fitness, heredity_fitness, a_Mp, b_Mp)

Heredity_To_F1 =  heredity_fitness < F1_fitness;

F2_offspring1 = F1_offspring(Heredity_To_F1,:);

r1= false(size(F2_offspring1));

nonSimilarBits = F2_offspring1 ~= heredity(Heredity_To_F1,:);

r1(nonSimilarBits) = rand(1, sum(sum(nonSimilarBits,2))) < a_Mp + (b_Mp -a_Mp)*rand(1, sum(sum(nonSimilarBits,2)));

F2_offspring1(r1) = heredity(r1);

F1_To_Heredity = F1_fitness < heredity_fitness ;

F2_offspring2 = heredity(F1_To_Heredity,:);

r2= false(size(F2_offspring2));

nonSimilarBits = F2_offspring2 ~= F1_offspring(F1_To_Heredity,:);

r2(nonSimilarBits) = rand(1, sum(sum(nonSimilarBits,2))) < a_Mp + (b_Mp -a_Mp)*rand(1, sum(sum(nonSimilarBits,2)));

F2_offspring2(r2) = F1_offspring(r2);

F1_offspring(Heredity_To_F1,:) = F2_offspring1;
F1_offspring(F1_To_Heredity,:) = F2_offspring2;

r(Heredity_To_F1,:)=r1;
r(F1_To_Heredity,:)=r2;

 Index1 = sum(r,2)>0;
% 
 feval = sum(Index1);

 F1_offspring = F1_offspring(Index1,:);


end
%% AntiSense Strand Generation
function[V, tempVf, tempVX, feval] = ASstrand(Id, V, OffspringSS, nVar, bits,...
            CostFunction, VarMin, VarMax, rr, selectVariable)
for j = 1: Id
   
  temp = sort(randperm(nVar, ceil(nVar*selectVariable))); 

  for i = temp
      if bits(i) ==3 % 3 is for binary variables
            k = 1;
            Offspring =OffspringSS(j, Bits);
                if rand > 0.5
                    V(j,sum(bits(1:i-1))+k)= 1-Offspring(k);
                else
                     V(j,sum(bits(1:i-1))+k)= Offspring(k);   
                end
      else
          % This calculation is for real and integer variables
      k = sort(randperm(bits(i), ceil(bits(i)*rr)));
             if i ==1
                Bits = 1:bits(1);    
             else
                Bits = sum(bits(1:i-1))+1:sum(bits(1:i));
            end
%       bits  = ((i - 1) * bits(i)) + 1 : i*bits(i);
      Offspring =fliplr(OffspringSS(j, Bits));
      V(j,sum(bits(1:i-1))+k)= Offspring(k);
      end
  end
 feval(j,1) = ~isequal(V(j,:),OffspringSS(j,:));
 end
feval = sum(feval);
[tempVf, tempVX] ...
         = variable_value(Id, V, bits, nVar, CostFunction, VarMin, VarMax);

end

%% Epimutation Process
function[heredity, feval] = Epimutation(v, b, mut_G_a, mut_G_b , id,...
                            heredity,  CostFunction,VarMin, VarMax, tau)
feval = 0;
counter = 0;
if ~isempty(id)
zeta = 4; % it is the research topic
while ~isempty(id) && (counter < zeta)
[affected, mutatedHeredity] = natural_mut(heredity.string(id,:),mut_G_a, mut_G_b, v, b, tau); 
 [Mut_f,  Mut_X] = variable_value(size(affected,1), affected, b, v, CostFunction, VarMin, VarMax); % Fitness and point calculation
 feval = feval + sum(mutatedHeredity);

 Id_best = heredity.f(id) > Mut_f;
updated = id(Id_best);
if isempty(updated)
counter = counter +1;
else
  if sum(Id_best) > 0   
      heredity.string(updated,:)= affected(Id_best,:);
      heredity.f(updated,:) = Mut_f (Id_best,:);
     heredity.X(updated,:)= Mut_X(Id_best,:);
  end  
end
id = setdiff(id, updated); % remove evolved candidates
end 
end
end
%% Mutation in Epimuation
 function [pop,  mutatedHeredity] = natural_mut(popu, mut_a, mut_b, nVar, bits, tau)
 pop = popu;
[s1, ~] = size(popu);
for i = 1: s1
    u = sort(randperm(nVar,ceil(tau*nVar)));
   for j = u
       	if j ==1
           	Bits = 1:bits(1);    
        else
            Bits = sum(bits(1:j-1))+1:sum(bits(1:j));
        end
       B = rand(1,bits(j)) < mut_a(j) + (mut_b(j) - mut_a(j))*rand(1,bits(j));
       if length(B)~= length(popu(i, Bits))
           a
       end
       pop(i, Bits) = xor(popu(i, Bits), B);
   end
    mutatedHeredity(i,1) = ~isequal(pop(i,:), popu(i,:)); 
end
 end