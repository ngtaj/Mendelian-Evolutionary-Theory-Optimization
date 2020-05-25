


%   23.	Non-Continuous Rotated Hybrid Composition Function 3
function fit=hybrid_rot_func3_noncont(x)
global aaa
persistent  o 
[ps,D]=size(x);
if aaa==0
    load hybrid_func3_data % saved the predefined optima, a 10*1000 matrix;
    o=o(1,1:D);
end
o=repmat(o,ps,1);
x=(abs(x-o)<0.5).*x+(abs(x-o)>=0.5).*(round(x.*2)./2);
fit=hybrid_rot_func3(x);
aaa=0;


function fit=hybrid_rot_func3(x)
global aaa
persistent  fun_num func o sigma lamda bias M
if aaa==0
    [ps,D]=size(x);
    aaa=1;
    fun_num=10;
    load hybrid_func3_data % saved the predefined optima, a 10*1000 matrix;
    if length(o(1,:))>=D
         o=o(:,1:D);
    else
         o=-5+10*rand(fun_num,D);
    end
    func.f1=str2func('fE_ScafferF6');
    func.f2=str2func('fE_ScafferF6');
    func.f3=str2func('frastrigin');
    func.f4=str2func('frastrigin');
    func.f5=str2func('fEF8F2');
    func.f6=str2func('fEF8F2');
    func.f7=str2func('fweierstrass');
    func.f8=str2func('fweierstrass');
    func.f9=str2func('fgriewank');
    func.f10=str2func('fgriewank');
    bias=((1:fun_num)-1).*100;
    sigma=[1,1,1,1,1,2,2,2,2,2];
    lamda=[5*5/100; 5/100; 5*1; 1; 5*1; 1; 5*10; 10; 5*5/200; 5/200];
    lamda=repmat(lamda,1,D);
    c=ones(1,D);
    if D==2,load hybrid_func3_M_D2,
    elseif D==10,load hybrid_func3_M_D10,
    elseif D==30,load hybrid_func3_M_D30,
    elseif D==50,load hybrid_func3_M_D50,
    else 
        for i=1:fun_num
           eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
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
    oo=repmat(o(i,:),ps,1);
    eval(['f=feval(func.f' int2str(i) ',((x-oo)./repmat(lamda(i,:),ps,1))*M.M' int2str(i) ');']);
    x1=5*ones(1,D);
    eval(['f1=feval(func.f' int2str(i) ',(x1./lamda(i,:))*M.M' int2str(i) ');']);
    fit1=2000.*f./f1;
    fit=fit+weight(:,i).*(fit1+bias(i));
end
%-------------------------------------------------
%basic functions

function f=fsphere(x)
%Please notice there is no use to rotate a sphere function, with rotation
%here just for a similar structure as other functions and easy programming
[ps,D]=size(x);
f=sum(x.^2,2);
%--------------------------------
function f=fsphere_noise(x)
[ps,D]=size(x);
f=sum(x.^2,2).*(1+0.1.*normrnd(0,1,ps,1));
%--------------------------------
function f=fgriewank(x)
[ps,D]=size(x);
f=1;
for i=1:D
    f=f.*cos(x(:,i)./sqrt(i));
end
f=sum(x.^2,2)./4000-f+1;
%--------------------------------
function f=fackley(x)
[ps,D]=size(x);
f=sum(x.^2,2);
f=20-20.*exp(-0.2.*sqrt(f./D))-exp(sum(cos(2.*pi.*x),2)./D)+exp(1);
%--------------------------------
function f=frastrigin(x)
[ps,D]=size(x);
f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
%--------------------------------
function f=frastrigin_noncont(x)
[ps,D]=size(x);
x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2)./2);
f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
%--------------------------------
function [f]=fweierstrass(x)
[ps,D]=size(x);
x=x+0.5;
a = 0.5;
b = 3;
kmax = 20;
c1(1:kmax+1) = a.^(0:kmax);
c2(1:kmax+1) = 2*pi*b.^(0:kmax);
f=0;
c=-w(0.5,c1,c2);
for i=1:D
f=f+w(x(:,i)',c1,c2);
end
f=f+c*D;

function y = w(x,c1,c2)
y = zeros(length(x),1);
for k = 1:length(x)
	y(k) = sum(c1 .* cos(c2.*x(:,k)));
end
%--------------------------------
function f=fE_ScafferF6(x)
fhd=str2func('ScafferF6');
[ps,D]=size(x);

f=0;
for i=1:(D-1)
    f=f+feval(fhd,(x(:,i:i+1)));
end
    f=f+feval(fhd,x(:,[D,1]));
%--------------------------------    
function f=fE_ScafferF6_noncont(x)
fhd=str2func('ScafferF6');
[ps,D]=size(x);
x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2)./2);
f=0;
for i=1:(D-1)
    f=f+feval(fhd,(x(:,i:i+1)));
end
    f=f+feval(fhd,x(:,[D,1]));
%------------------------------
function f=fEF8F2(x)
[ps,D]=size(x);
f=0;
for i=1:(D-1)
    f=f+F8F2(x(:,[i,i+1]));
end
    f=f+F8F2(x(:,[D,1]));

%--------------------------------
function f=fschwefel_102(x)
[ps,D]=size(x);
f=0;
for i=1:D
    f=f+sum(x(:,1:i),2).^2;
end
%--------------------------------
function f=felliptic(x)
[ps,D]=size(x);
a=1e+6;
f=0;
for i=1:D
f=f+a.^((i-1)/(D-1)).*x(:,i).^2;
end
function f=ScafferF6(x)
f=0.5+(sin(sqrt(x(:,1).^2+x(:,2).^2)).^2-0.5)./(1+0.001*(x(:,1).^2+x(:,2).^2)).^2;
function f=F8F2(x)
f2=100.*(x(:,1).^2-x(:,2)).^2+(1-x(:,1)).^2;
f=1+f2.^2./4000-cos(f2);
  


