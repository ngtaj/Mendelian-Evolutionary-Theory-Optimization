%% Functions to test
function [CostFunction,  VarMin, VarMax] = TestFunctionsMETOPaper(func_num, D)
 fun_num =10; % ten functions are taken for composition 
%% Multimodal functions
if func_num == 1; CostFunction = @(x)  InvertedCosineWave(x);   VarMin = -5;   VarMax = 5;   end
if func_num == 2;  CostFunction = @(x) rastr(x);   VarMin=-5.125;  VarMax= 5.125;  end 
if func_num == 3; CostFunction = @(x)  GenSchwefelNo226(x);  VarMin = -500;   VarMax = 500;   end 
if func_num == 4;  CostFunction = @(x) wavy(x);   VarMin = -pi;  VarMax = pi;  end
if func_num == 5;  CostFunction = @(x) stochasticf(x);  VarMin = -5;   VarMax = 5;  end
if func_num == 6;  CostFunction = @(x) dropwave(x);   VarMin = -5.125;  VarMax = 5.125;  end
if func_num == 7;  CostFunction = @(x) LangMann(x); VarMin = 0;  VarMax = 10;   end
if func_num == 8; CostFunction = @(x)  DeflectedCorrugatedSpring(x);   VarMin = 0;   VarMax = 10;   end
if func_num == 9;  CostFunction = @(x) LunacekBiRastrigin(x); VarMin =-5.125;  VarMax = 5.125;   end 
if func_num == 10;  CostFunction = @(x) shubert3fcn(x);   VarMin = -10;  VarMax =  10; end
if func_num == 11;  CostFunction = @(x) shubert4fcn(x);   VarMin = -10;  VarMax =  10; end
if func_num == 12;  CostFunction = @(x) crossintray(x);   VarMin = -10;  VarMax =  10; end
if func_num == 13;  CostFunction = @(x) holdertable(x);   VarMin = -10;  VarMax =  10; end
if func_num == 14;  CostFunction = @(x) bird(x);   VarMin = -10;  VarMax =  10; end 
if func_num == 15;  CostFunction = @(x) periodicfcn(x);   VarMin = -5;  VarMax =  10;  end
if func_num == 16;  CostFunction = @(x) Deb01(x);   VarMin = -1;  VarMax =  1; end 
if func_num == 17;  CostFunction = @(x) eggholder(x);  VarMin = -512;     VarMax = 512;  end  
if func_num == 18;CostFunction = @(x) deceptive (x);  VarMin=0;  VarMax= 1; end
if func_num == 19;  CostFunction = @(x) michal_Ext(x); VarMin = 0;   VarMax = pi;  end      
if func_num == 20; CostFunction = @(x) LorF2(x); VarMin = 0;   VarMax = 1;  end 
if func_num == 21;  CostFunction = @(x) deformedschaffer2(x);   VarMin = -100;  VarMax = 100;  end
if func_num == 22; CostFunction = @(x) gramacyleefcn(x); VarMin =-0.5; VarMax = 2.5; end
if func_num == 23;   CostFunction = @(x) xinsheyangn4fcn(x);   VarMin = -100;  VarMax =  100; end
if func_num == 24;  CostFunction = @(x) AlpineN01(x);   VarMin = -10;  VarMax =  10; end
if func_num == 25;  CostFunction = @(x) weierstrassfunc(x);   VarMin = -0.5;  VarMax = 0.5;  end
if func_num == 26;  CostFunction = @(x) giunta(x);  VarMin = -1;   VarMax = 1;  end
if func_num == 27;  CostFunction = @(x) mixSimpletonSquare(x);  VarMin = 0;   VarMax = 10;  end
if func_num == 28;  CostFunction = @(x) sumOfSquare(x);  VarMin = 0;   VarMax = 10;  end
if func_num == 29;  CostFunction = @(x) Noisy_SinusoidalFunction(x);  VarMin = -100;   VarMax = 100;  end
if func_num == 30;  CostFunction = @(x) schaffer2(x); VarMin = -100; VarMax = 100;  end
if func_num == 31;  CostFunction = @(x) schaffer4(x); VarMin = -100; VarMax = 100;  end
if func_num == 32;  CostFunction = @(x) schaffer6(x); VarMin = -100; VarMax = 100;  end
if func_num == 33;  CostFunction = @(x) rosenbrock_Leon(x); VarMin = -1.2; VarMax = 1.2;  end
if func_num == 34; CostFunction = @(x)  DixonPrice (x);  VarMin = -10;     VarMax = 10;  end
if func_num == 35;  CostFunction = @(x) h11(x);  VarMin = 3;   VarMax = 13;  end
if func_num == 36;  CostFunction = @(x) griewank(x); VarMin = -100;  VarMax = 100;  end
if func_num == 37;  CostFunction = @(x) weierstrassfunc(x);   VarMin = -0.5;  VarMax = 0.5;  end

%% Discontinuous functions:

if func_num == 38;  CostFunction = @(x) frastrigin_noncont(x);   VarMin=-5.125;  VarMax= 5.125;  end  
if func_num == 39;  CostFunction = @(x)   fE_ScafferF6_noncont(x); VarMin = -100; VarMax = 100;  end 

%% Multi-objective functions:
 if func_num == 40;  CostFunction = @(x)   kursaweFcn(x);   VarMin = -5;  VarMax = 5;  end

%% Noisy functions:
 
if func_num == 41;  CostFunction = @(x) Noisy_styblinskitang(x);  VarMin = -5;   VarMax = 5;  end
if func_num == 42;  CostFunction = @(x) eggholderNoisy(x);   VarMin = -512;  VarMax =  512; end 
if func_num == 43;  CostFunction = @(x)  LangMannNoisy(x); VarMin = 0;  VarMax = 10;   end
if func_num == 44;  CostFunction = @(x) Noisy_deceptive(x);  VarMin = 0;   VarMax = 1;  end

%% Constrained functions: 
if func_num == 45;  CostFunction = @(x) PrG3f(x);  VarMin = 0;   VarMax = 10;  end
if func_num == 46; CostFunction = @(x) keane94Bump(x); VarMin=0; VarMax= 10; end
% 
% %% Integrated functions:
if func_num == 47;  CostFunction = @(x)  mixSimpletonSquare(x);  VarMin = 0;   VarMax = 10;  end

if func_num == 48;  CostFunction = @(x) IntegratedFunc4(x);   VarMin = -pi;  VarMax = pi;  end
if func_num == 49;  CostFunction = @(x) IntegratedFunc5(x);   VarMin = -10;  VarMax = 10;  end
if func_num == 50; CostFunction = @(x) IntegratedFunc11(x);   VarMin = -10;  VarMax = 10;  end
if func_num == 51; CostFunction = @(x) IntegratedFunc10(x);   VarMin = -5;  VarMax = 5;  end

%% Rotated Functions

if func_num == 52; VarMin =-5; VarMax = 5; load rastrigin_func_data;
        CostFunction = @(x) rastrigin_func(x, o, M, A, B, Sigma, lamda, bias, fun_num ); end


%% Hybrid Composite funciton 1

o =[]; M=[]; A=[]; B=[]; sigma=[]; lamda=[]; bias=[];

if func_num == 53
    VarMin =-5; VarMax = 5; 
    load hybrid_func1_data;
        bias=((1:fun_num)-1).*100; 
        Sigma=ones(1,fun_num);
        lamda=[1; 1; 10; 10; 5/60; 5/60; 5/32; 5/32; 5/100; 5/100];
        lamda=repmat(lamda,1,D);
           for i=1:fun_num
                 eval(['M.M' int2str(i) '=diag(ones(1,D));']);
           end
    CostFunction = @(x) hybrid_func1(x, o, M, A, B, Sigma, lamda, bias, fun_num ); 
end
%% Hybrid Composite funciton 2
if func_num == 54  
    VarMin =-5; VarMax = 5; 
            load hybrid_func1_data;
            bias=((1:fun_num)-1).*100;
            Sigma=ones(1,fun_num); 
            lamda=[1; 1; 10; 10; 5/60; 5/60; 5/32; 5/32; 5/100; 5/100];
            lamda=repmat(lamda,1,D);
            c=[2,2,2,2,2,2,2,2,2,2,2];
            load hybrid_func1_M_D50; % To test 30 variables problem change data file to "hybrid_func1_M_D30" 
    CostFunction = @(x) hybrid_rot_func1(x, o, M, A, B, Sigma, lamda, bias, fun_num ); 
end
%% Hybrid Composite funciton 3
if func_num == 55 
    VarMin =-5; VarMax = 5; 
            load hybrid_func1_data;
            bias=((1:fun_num)-1).*100;
            Sigma=ones(1,fun_num); 
            lamda=[1; 1; 10; 10; 5/60; 5/60; 5/32; 5/32; 5/100; 5/100];
            lamda=repmat(lamda,1,D);
            c=[2,2,2,2,2,2,2,2,2,2,2];
            load hybrid_func1_M_D50; % To test 30 variables problem change data file to "hybrid_func1_M_D30" 

    
    CostFunction = @(x) hybrid_rot_func1_noise(x, o, M, A, B, Sigma, lamda, bias, fun_num ); 
end

%% Hybrid Composite funciton 4
if func_num == 56 
    VarMin =-5; VarMax = 5;   
            load hybrid_func2_data; 
            bias=((1:fun_num)-1).*100;
            Sigma=[1 2 1.5 1.5 1 1 1.5 1.5 2 2];
            lamda=[2*5/32; 5/32; 2*1; 1; 2*5/100; 5/100; 2*10; 10; 2*5/60; 5/60];
            lamda=repmat(lamda,1,D);
            c=[2 3 2 3 2 3 20 30 200 300];
            load hybrid_func2_M_D50; % To test 30 variables problem change data file to "hybrid_func1_M_D30" 

    CostFunction = @(x) hybrid_rot_func2(x, o, M, A, B, Sigma, lamda, bias, fun_num ); 
end

%% Hybrid Composite funciton 5    

if func_num == 57
    VarMin =-5; VarMax = 5 ;   
        load hybrid_func2_data;  
        
        bias=((1:fun_num)-1).*100;
    Sigma=[0.1 2 1.5 1.5 1 1 1.5 1.5 2 2];
    lamda=[0.1*5/32; 5/32; 2*1; 1; 2*5/100; 5/100; 2*10; 10; 2*5/60; 5/60];
    lamda=repmat(lamda,1,D);
    c=[2 3 2 3 2 3 20 30 200 300];
load hybrid_func2_M_D50; % To test 30 variables problem change data file to "hybrid_func1_M_D30" 

    CostFunction = @(x) hybrid_rot_func2_narrow(x, o, M, A, B, Sigma, lamda, bias, fun_num ); 
end

%% Hybrid Composite funciton 6
if func_num == 58 
    VarMin =-5; VarMax = 5; 
            load hybrid_func2_data;
            bias=((1:fun_num)-1).*100;
            Sigma=[1 2 1.5 1.5 1 1 1.5 1.5 2 2];
            lamda=[2*5/32; 5/32; 2*1; 1; 2*5/100; 5/100; 2*10; 10; 2*5/60; 5/60];
             lamda=repmat(lamda,1,D);
            load hybrid_func2_M_D50;% To test 30 variables problem change data file to "hybrid_func1_M_D30" 

    CostFunction = @(x) hybrid_rot_func2_onbound(x, o, M, A, B, Sigma, lamda, bias, fun_num ); 
end


%% Hybrid Composite funciton 7
if func_num == 59 
    VarMin =-5; VarMax = 5; 
        load hybrid_func3_data ; 
        
        bias=((1:fun_num)-1).*100;
    Sigma=[1,1,1,1,1,2,2,2,2,2];
    lamda=[5*5/100; 5/100; 5*1; 1; 5*1; 1; 5*10; 10; 5*5/200; 5/200];
    lamda=repmat(lamda,1,D);
    c=ones(1,D);
    load hybrid_func3_M_D50;% To test 30 variables problem change data file to "hybrid_func1_M_D30" 

    CostFunction = @(x) hybrid_rot_func3(x, o, M, A, B, Sigma, lamda, bias, fun_num ); 
end

end
%% Functions to use in Composite Funcitons

% 
function f=rastrigin_func(x, o, M, A, B, sigma, lamda, bias, fun_num)
[ps,D]=size(x);
o = o(1:D);
x=x-repmat(o,ps,1);
f  = feval(@rastr,x);
% f = awgn(f,10);
end
% 

function fit=hybrid_func1(x, o, M, A, B, sigma, lamda, bias, fun_num)
persistent func
    [ps,D]=size(x);
      o = o(:,1:D);
    func.f1=str2func('rastr');
    func.f2=str2func('rastr');
    func.f3=str2func('weierstrassfunc');
    func.f4=str2func('weierstrassfunc');
    func.f5=str2func('griewank');
    func.f6=str2func('griewank');
    func.f7=str2func('ackley');
    func.f8=str2func('ackley');
    func.f9=str2func('spheref');
    func.f10=str2func('spheref');
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
end


function fit=hybrid_rot_func1(x, o, M, A, B, sigma, lamda, bias, fun_num)
persistent  func 
[ps,D]=size(x);
      o = o(:,1:D);
   func.f1=str2func('rastr');
    func.f2=str2func('rastr');
    func.f3=str2func('weierstrassfunc');
    func.f4=str2func('weierstrassfunc');
    func.f5=str2func('griewank');
    func.f6=str2func('griewank');
    func.f7=str2func('ackley');
    func.f8=str2func('ackley');
    func.f9=str2func('spheref');
    func.f10=str2func('spheref');
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
end

function fit=hybrid_rot_func2(x, o, M, A, B, sigma, lamda, bias, fun_num)
persistent  func
[ps,D]=size(x);
      o = o(:,1:D);
    o(10,:)=0;
    func.f1=str2func('ackley');
    func.f2=str2func('ackley');
    func.f3=str2func('rastr');
    func.f4=str2func('rastr');
    func.f5=str2func('spheref');
    func.f6=str2func('spheref');
    func.f7=str2func('weierstrassfunc');
    func.f8=str2func('weierstrassfunc');
    func.f9=str2func('griewank');
    func.f10=str2func('griewank');
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
end

function fit=hybrid_rot_func2_narrow(x, o, M, A, B, sigma, lamda, bias, fun_num)
persistent  func
[ps,D]=size(x);
      o = o(:,1:D);
    o(10,:)=0;
     func.f1=str2func('ackley');
    func.f2=str2func('ackley');
    func.f3=str2func('rastr');
    func.f4=str2func('rastr');
    func.f5=str2func('spheref');
    func.f6=str2func('spheref');
    func.f7=str2func('weierstrassfunc');
    func.f8=str2func('weierstrassfunc');
    func.f9=str2func('griewank');
    func.f10=str2func('griewank');
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
end

function fit=hybrid_rot_func2_onbound(x, o, M, A, B, sigma, lamda, bias, fun_num)
persistent func
[ps,D]=size(x);
      o = o(:,1:D);
    o(10,:)=0;
    o(1,2.*[1:floor(D/2)])=5;
    func.f1=str2func('ackley');
    func.f2=str2func('ackley');
    func.f3=str2func('rastr');
    func.f4=str2func('rastr');
    func.f5=str2func('spheref');
    func.f6=str2func('spheref');
    func.f7=str2func('weierstrassfunc');
    func.f8=str2func('weierstrassfunc');
    func.f9=str2func('griewank');
    func.f10=str2func('griewank');
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
end

function fit=hybrid_rot_func3(x, o, M, A, B, sigma, lamda, bias, fun_num)
persistent func 
[ps,D]=size(x);
      o = o(:,1:D);
  o=-5+10*rand(fun_num,D);
    func.f1=str2func('fE_ScafferF6');
    func.f2=str2func('fE_ScafferF6');
    func.f3=str2func('rastr');
    func.f4=str2func('rastr');
    func.f5=str2func('fEF8F2');
    func.f6=str2func('fEF8F2');
    func.f7=str2func('weierstrassfunc');
    func.f8=str2func('weierstrassfunc');
    func.f9=str2func('griewank');
    func.f10=str2func('griewank');
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
end


function fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M)
[ps,D]=size(x);
for i=1:fun_num
    oo=repmat(o(i,:),ps,1);
    weight(:,i)=exp(-sum((x-oo).^2,2)./2./(D*sigma(i)^2));
end

[tmp,tmpid]=sort(weight,2);
for i=1:ps
    weight(i,:)=(weight(i,:)==tmp(i,fun_num)).*weight(i,:)+(weight(i,:)~=tmp(i,fun_num)).*(weight(i,:).*(1-tmp(i,fun_num).^10));
end
weight=weight./repmat(sum(weight,2),1,fun_num);

fit=0;
for i=1:fun_num
%     oo=repmat(o(i,:),ps,1);
    eval(['f=feval(func.f' int2str(i) ',((x-oo)./repmat(lamda(i,:),ps,1))*M.M' int2str(i) ');']);
    x1=5*ones(1,D);
    eval(['f1=feval(func.f' int2str(i) ',(x1./lamda(i,:))*M.M' int2str(i) ');']);
    fit1=2000.*f./f1;
    fit=fit+weight(:,i).*(fit1+bias(i));
end
end

%% Basic functions

function f=fE_ScafferF6_noncont(x)
fhd=str2func('schaffer6');
[ps,D]=size(x);
x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2)./2);
f=0;
for i=1:(D-1)
    f=f+feval(fhd,(x(:,i:i+1)));
end
    f=f+feval(fhd,x(:,[D,1]));
end
%------------------------------
function f=fEF8F2(x)
[ps,D]=size(x);
f=0;
for i=1:(D-1)
    f=f+F8F2(x(:,[i,i+1]));
end
    f=f+F8F2(x(:,[D,1]));
end

function f=F8F2(x)
f2=100.*(x(:,1).^2-x(:,2)).^2+(1-x(:,1)).^2;
f=1+f2.^2./4000-cos(f2);
end
function f=frastrigin_noncont(x)
[ps,D]=size(x);
x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2)./2);
f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
end
function f=fsphere_noise(x)
[ps,D]=size(x);
f=sum(x.^2,2).*(1+0.1.*randn(ps,1));
end

function f=felliptic(x)
[ps,D]=size(x);
a=1e+6;
f=0;
for i=1:D
f=f+a.^((i-1)/(D-1)).*x(:,i).^2;
end
end
function f=fE_ScafferF6(x)
fhd=str2func('schaffer6');
[ps,D]=size(x);

f=0;
for i=1:(D-1)
    f=f+feval(fhd,(x(:,i:i+1)));
end
    f=f+feval(fhd,x(:,[D,1]));
end
