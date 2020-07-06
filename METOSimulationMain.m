clc
clearvars

%  In the following addpath all files
% and Testfunctions folder in in E drive. In your system drive letter will
% chage. Please chnage the path of METO Files before start optimizing the problem.
addpath 'E:\METOforGitHub'
addpath 'E:\METOforGitHub\Testfunctions'
global vType
Iteration = 250000; % Number of evolution epochs 
max_f_evaluation =100000;

%% Test problems:
nPop= 100; % number of individuals in a populaiton
nVar =100; % Number of variables in the optimizaiton Problem


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
func_num = 3; % Change the Function number of 1 to 53
% problems 54 to 59 are for 50 variables only. For these problems you have to assign "nVar =50", above. 
% These functions can run on 30 variables, If you replace data files from "hybrid_*****_M_D50" to "hybrid_*****_M_D30" in the "TestFunctionsMETOPaper.m" file. 
%These data files are taken from CEC2014 and can be see in "Testfunctions" folder. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Select Benchmark Function for optimizing
[CostFunction,  VarMin, VarMax] = TestFunctionsMETOPaper(func_num, nVar);
% Run METO Algorithm by following command

% Defining the lower and upper bound for each variables
if length(VarMin) < 2 % Upper and Lower Bound on nVar variables
   VarMin = repmat(VarMin,1,nVar); 
   VarMax = repmat(VarMax,1,nVar); 
end


%% Dealing with the type of Variables
vType = ones(1,nVar); % In this vector, replace 1 by either 2 or 3 according to the problem
% vType is a vector of [1 1 1 1 ... 1], all "1"  indicates that all variables are Real 
% 1 = Real variables ==> (vType = [1 1 1 1 ... 1])
% 2 = Integer variables ==>  (vType = [1 2 2 1 2... 2], here  2, 3, 5
% and nth variables are Integer type)
% 3 = Binary variables ==>(vType = [1 3 2 1 3... 2], here  2, 5 variables are
% Binary type)
% If all variables are Binary ==>(vType = [3 3 3 3 ... 3])
% If all variables are Integer ==>(vType = [2 2 2 2 ... 2])

%% Run METO
   [BestCost, ConvergenceCurve, FunctionEvaluation] = METO(CostFunction,nVar, VarMin, VarMax, vType, nPop, Iteration, max_f_evaluation);
 figure
 plot(ConvergenceCurve,'LineWidth',2);
 ylabel('Function Value');