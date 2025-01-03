clearvars
clc
load("colour.mat")
format long 
syms x 
delta =1;
rho=-2.8;
a2=-1;
eps= 1;
lambda_s=-5.9;
lambda_u=4.5;
a1=-1;
mu1=-0.1;
mu2=0.1;
%% Write T12
% mu1 = 0
if mu1==0
    A0=delta*exp(-rho/delta);
   T12=@(x) A0*exp(rho/x);
% mu1 > 0
elseif mu1>0
    A1=delta*exp((rho/(2*sqrt(mu1)))*(log(abs(delta/sqrt(mu1)-1))-log(abs(delta/sqrt(mu1)+1))));
   T12=@(x) A1*exp((rho/(2*sqrt(mu1)))*(-log(abs(x/sqrt(mu1)-1))+log(abs(x/sqrt(mu1)+1))));
% mu > 0
else   
   A2 = delta*exp((rho/sqrt(abs(mu1)))*atan(delta/sqrt(abs(mu1))));
  T12=@(x)A2* exp(-(rho/sqrt(abs(mu1)))*atan(x/sqrt(abs(mu1))));
end
%% Write T34
v = lambda_s/lambda_u; 
D = -eps*abs(-eps)^v;
T34 = @(x) D * abs(x)^(-v);
%%
close all
f1=figure(2);
f1.Units="centimeters";
f1.OuterPosition = [2 10 23 17];
hold on; box on; grid off
set ( gca , 'FontSize' , 12 , 'fontname' , 'times');
xlim([-1 1])
ylim([-1 1])
xlabel('$\mu_2$', Interpreter='latex')
ylabel('$\mu_3$', Interpreter='latex')
%% Homoclinic
y3 = linspace(-1,1,100000);
for i =1:length(y3)
    x3(i)=-a1*T12(y3(i));
end
plot(x3,y3,"LineWidth",1.5,'Color',colour.grass)
legend('HC $p_2 \rightarrow p_2$','FontSize' , 11, 'interpreter', 'latex',    'Location','northeastoutside','FontSize' , 11, 'interpreter', 'latex')
% saveas(f1, 'het.svg', 'svg')



















