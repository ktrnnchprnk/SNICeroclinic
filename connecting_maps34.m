clearvars; clc
close all
load("colour.mat")
format long 
syms x
eps= 1;
lambda_s=-4;
lambda_u=2;
delta =1;
rho=-2.8;
mu3=0;
v = lambda_s/lambda_u; 
D = -eps*abs(-eps)^v;
T23 = @(x) D * abs(x)^(-v);
n=10;
x0 = linspace(-n,0,100000);
for i =1:length(x0)
    y1(i)=T23(x0(i));
end
A0=delta*exp(-rho/delta);
T41=@(x) A0*exp(rho/x);
f1=figure(1);
f1.Units="centimeters";
f1.OuterPosition = [2 10 23 10];
subplot(1,3,1)
hold on; box on; grid on
set ( gca , 'FontSize' , 12 , 'FontName', 'Times New Roman');
subtitle('$|\lambda_s|>|\lambda_u|$',Interpreter='latex')
ylabel('$T_{34}(x)$', Interpreter='latex')
xlabel('$x$', Interpreter='latex')
title('A')
scatter(x0,y1,5,'filled','o','k')
axis on
ylim([-2 1])
xlim([-2 1])
%%
eps= 1;
lambda_s=-2;
lambda_u=2;
delta =1;
rho=-2.8;
mu3=0;
v = lambda_s/lambda_u; 
D = -eps*abs(-eps)^v;
T23 = @(x) D * abs(x)^(-v);
n=10;
x0 = linspace(-n,0,100000);
for i =1:length(x0)
    y1(i)=T23(x0(i));
end
subplot(1,3,2)
hold on; box on; grid on
subtitle('$|\lambda_s|=|\lambda_u|$',Interpreter='latex')
set ( gca , 'FontSize' , 12 , 'fontname' , 'times');
ylabel('$T_{34}(x)$', Interpreter='latex')
xlabel('$x$', Interpreter='latex')
scatter(x0,y1,5,'filled','o','k')
% plot([0, 0], [-1 10],LineWidth=1,Color=colour.grey)
title('B')
axis on
ylim([-2 1])
xlim([-2 1])
%%
eps= 1;
lambda_s=-2;
lambda_u=5;
delta =1;
rho=-2.8;
mu3=0;
v = lambda_s/lambda_u; 
D = -eps*abs(-eps)^v;
T23 = @(x) D * abs(x)^(-v);
n=10;
x0 = linspace(-n,0,100000);
for i =1:length(x0)
    y1(i)=T23(x0(i));
end
subplot(1,3,3)
hold on; box on; grid on
subtitle('$|\lambda_s|<|\lambda_u|$',Interpreter='latex')
set ( gca , 'FontSize' , 12 , 'fontname' , 'times');
ylabel('$T_{34}(x)$', Interpreter='latex')
xlabel('$x$', Interpreter='latex')
title('C')
scatter(x0,y1,5,'filled','o','k')
axis on
ylim([-2 1])
xlim([-2 1])
% saveas(f1, 'connecting_maps34.svg', 'svg')
% print(f1, 'connecting_maps34.eps', '-depsc')