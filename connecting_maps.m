clearvars; clc
close all

addpath("/home/kn356/Desktop/FILES/StandardMATLAB/")
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

n=4;
x0 = linspace(-n,0,100000);
for i =1:length(x0)
    y1(i)=T23(x0(i));
end

A0=delta*exp(-rho/delta);
T41=@(x) A0*exp(rho/x);

f1=figure(1);
f1.Units="centimeters";
f1.OuterPosition = [2 10 23 9];


n=5;
x0 = linspace(0,n,100000);
for i =1:length(x0)
    y2(i)=T41(x0(i));
end

subplot(1,3,3)
hold on; box on; grid on
set ( gca , 'FontSize' , 12 , 'fontname' , 'times');
ylabel('$T_{12}(x)$', Interpreter='latex')
xlabel('$x$', Interpreter='latex')
% fill([-1.96 -1.96 0 0 -1.96], [-0.9 9.93 9.93 -0.93 -0.93], colour.grey, 'EdgeColor', 'none' )
% plot([-1000, 1000], [0 0],LineWidth=1,Color='k',LineStyle=':')
% plot([0 0],[-1000, 1000],LineWidth=1,Color='k',LineStyle=':')
scatter(x0,y2,5,'filled','o','k')

% plot([0, 0], [-1 10],LineWidth=1,Color=colour.grey)
subtitle('$\mu_1=0$',Interpreter='latex')
axis on
ylim([-1 10])
xlim([-2 n])
title('C')

mu3=0.1;
A1=delta*exp((rho/(2*sqrt(mu3)))*(log(abs(delta/sqrt(mu3)-1))-log(abs(delta/sqrt(mu3)+1))));
T41=@(x) A1*exp((rho/(2*sqrt(mu3)))*(-log(abs(x/sqrt(mu3)-1))+log(abs(x/sqrt(mu3)+1))));
n=5;
x0 = linspace(sqrt(mu3),n,100000);
for i =1:length(x0)
    y3(i)=T41(x0(i));
end

%%
subplot(1,3,2)
hold on; box on; grid on
set ( gca , 'FontSize' , 12 , 'fontname' , 'times');
ylabel('$T_{12}(x)$', Interpreter='latex')
xlabel('$x$', Interpreter='latex')

% fill([-1.96 -1.96 0 0 -1.96], [-0.9 9.93 9.93 -0.93 -0.93], colour.grey, 'EdgeColor', 'none' )
% plot([-1000, 1000], [0 0],LineWidth=1,Color='k',LineStyle=':')
% plot([0 0],[-1000, 1000],LineWidth=1,Color='k',LineStyle=':')
scatter(x0,y3,5,'filled','o','k')
plot([sqrt(mu3), sqrt(mu3)], [-10 10], 'LineWidth',2, 'LineStyle',':', 'Color',colour.grey)


subtitle('$\mu_1>0$',Interpreter='latex')
axis on
ylim([-1 10])
xlim([-2 n])
title('B')

mu3=-0.1;
A2 = delta*exp((rho/sqrt(abs(mu3)))*atan(delta/sqrt(abs(mu3))));
T41=@(x)A2* exp(-(rho/sqrt(abs(mu3)))*atan(x/sqrt(abs(mu3))));
n=5;
x0 = linspace(-n,n,100000);
for i =1:length(x0)
    y4(i)=T41(x0(i));
end


subplot(1,3,1)
hold on; box on; grid on
set ( gca , 'FontSize' , 12 , 'fontname' , 'times');
ylabel('$T_{12}(x)$', Interpreter='latex')
xlabel('$x$', Interpreter='latex')

% fill([-1.96 -1.96 0 0 -1.96], [-0.9 9.93 9.93 -0.93 -0.93], colour.grey, 'EdgeColor', 'none' )
% plot([-1000, 1000], [0 0],LineWidth=1,Color='k',LineStyle=':')
% plot([0 0],[-1000, 1000],LineWidth=1,Color='k',LineStyle=':')
scatter(x0,y4,5,'filled','o','k')
% plot(x0, x0,LineWidth=1,Color=colour.grey)
subtitle('$\mu_1<0$',Interpreter='latex')
title('A')
axis on
ylim([-1 10])
xlim([-2 n])
set(gcf, 'Renderer', 'Painters');
saveas(f1, 'connecting_maps.svg', 'svg')

print(f1, 'connecting_maps.eps', '-depsc')

print(f1, 'connecting_maps', '-dpng', '-r1200')








