clc; clearvars
load("colour.mat")
load('params.mat')
% Define parameters where het orbit is
p.beta=5.16184E-03;  p.b =2.53003E-01;
p.phi=0.9;p.phi2=0.8;
p.m=4;
% Compute equilibria and corresponding eigenvalues and -vectors 
[u,e,v1,v2] = compute_fp(@(x) GTP2([], x,p),[0.001, 1.2],[0.001, 1],0.02,0.02);
%% Compute the separatrices for fp1
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
epsi=[0.01,0.01];
X0s1 = u(3,:)+epsi.*v2(3,:);
X0s2 = u(3,:)-epsi.*v1(3,:);
X0u1 = u(3,:)+epsi.*v1(3,:);
X0u2 = u(3,:)+epsi.*v2(3,:);
Tspanb1 = 1000:-0.1:0;
Tspanb2 = 90:-0.1:0;
Tspanf = 0:0.1:100;
Tspanf2 = 0:0.1:20;
[T, Y1]=ode45(@GTP2,Tspanb1,X0s1,opts, p);
[~, Y2]=ode45(@GTP2,Tspanf,X0u1,opts, p);
[~, Y3]=ode45(@GTP2,Tspanb2,X0s2,opts, p);
[~, Y4]=ode45(@GTP2,Tspanf2,X0u2,opts, p);
%% Compute the separatrices for fp2
opts = odeset('RelTol',1e-4,'AbsTol',1e-6);
epsi=[0.001,0.005];
X0s1 = [0.1 0.92];
X0s2 = u(1,:)-epsi.*v1(1,:);
X0u1 =[0.225 0.932];
X0u2 = [0.3 0.3];
Tspanf = 0:0.01:1;
[~, Y12]=ode45(@GTP2,Tspanf,X0s1,opts, p);
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
Tspanf = 500:-0.1:0;
[T, Y22]=ode45(@GTP2,Tspanf,X0u1,opts, p);
%% Plot the phase portrait
close all
f2=figure(2);
f2.Units="centimeters";
f2.OuterPosition = [2 10 15 12];
hold on; box on; grid off
set ( gca , 'FontSize' , 12 , 'fontname' , 'times');
ylabel('$L$', Interpreter='latex')
xlabel('$G$', Interpreter='latex')
N=9000;
% Plot separatrices
plot(Y1(1:N,1),Y1(1:N,2), 'Color','k', 'LineWidth',1.5, 'LineStyle','--')
plot(Y2(:,1),Y2(:,2), 'Color',colour.grey, 'LineWidth',1.5)
plot(Y3(:,1),Y3(:,2), 'Color',colour.grey, 'LineWidth',1.5)
plot(Y12(:,1),Y12(:,2), 'Color',colour.grey, 'LineWidth',1.5)
plot(Y22(:,1),Y22(:,2), 'Color',colour.grey, 'LineWidth',1.5)
% Plot fixed points
plot(u(1,1),u(1,2),'o','MarkerFaceColor', ...
colour.yellow,'MarkerEdgeColor','black', 'MarkerSize', 8,'LineWidth',1)
plot(u(3,1),u(3,2),'o','MarkerFaceColor', ...
colour.pink,'MarkerEdgeColor','black', 'MarkerSize', 8,'LineWidth',1)
xlim([0.15 0.7])
ylim([0.2 1])
% print(f2, 'GTP_pp.eps', '-depsc')
% saveas(f2,'GTP_pp.svg', 'svg')