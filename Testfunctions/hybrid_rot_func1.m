function fit=hybrid_rot_func1(x)


persistent  fun_num func o sigma lamda bias M
if aaa==0
    [ps,D]=size(x);
    aaa=1;
    fun_num=10;
    load hybrid_func1_data % saved the predefined optima
    if length(o(1,:))>=D
         o=o(:,1:D);
    else
         o=-5+10*rand(fun_num,D);
    end
    func.f1=str2func('frastrigin');
    func.f2=str2func('frastrigin');
    func.f3=str2func('fweierstrass');
    func.f4=str2func('fweierstrass');
    func.f5=str2func('fgriewank');
    func.f6=str2func('fgriewank');
    func.f7=str2func('fackley');
    func.f8=str2func('fackley');
    func.f9=str2func('fsphere');
    func.f10=str2func('fsphere');
    bias=((1:fun_num)-1).*100;
    sigma=ones(1,fun_num); 
    lamda=[1; 1; 10; 10; 5/60; 5/60; 5/32; 5/32; 5/100; 5/100];
    lamda=repmat(lamda,1,D);
    c=[2,2,2,2,2,2,2,2,2,2,2];
    if D==2,load hybrid_func1_M_D2,
    elseif D==10,load hybrid_func1_M_D10,
    elseif D==30,load hybrid_func1_M_D30,
    elseif D==50,load hybrid_func1_M_D50,
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
function [q, r] = cGram_Schmidt (A)
% computes the QR factorization of $A$ via
% classical Gram Schmid
%
[n, m] = size(A);
q = A;
for j = 1 : m
    for i = 1 : j - 1
        r(i, j) = q(:, j)' * q(:, i);
    end
    for i = 1 : j - 1
        q(:, j) = q(:, j) - r(i, j) * q(:, i);
    end
    t = norm(q(:, j), 2) ;
    q(:, j) = q(:, j) / t ;
    r(j, j) = t ;
end

function M = rot_matrix(D, c)
A = normrnd(0, 1, D, D);
P = cGram_Schmidt(A);
A = normrnd(0, 1, D, D);
Q = cGram_Schmidt(A);
u = rand(1, D);
D = c .^ ((u - min(u)) ./ (max(u) - min(u)));
D = diag(D);
M = P * D * Q;