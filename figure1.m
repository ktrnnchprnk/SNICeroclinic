clc; clearvars
% Load colour matrix
load("colour.mat")
% Set the parameters
p.a=0.042; p.b=0.49577; p.c=-0.85; p.eps=1;
% Compute fixed points and the corresponding eigenvalues and eigenvectors
[u,e,v1,v2] = compute_fp(@(x) SNICeroclinic([], x,p),[-6, 6],[-6 6],0.1,0.1);
%% Compute the separatrices #1
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
epsi=[0.1,0.1];
X0u2 = u(1,:)+epsi.*v1(1,:);
Tspanf2 = 0:0.1:1000;
[~, Y1]=ode15s(@SNICeroclinic,Tspanf2,X0u2,opts, p);
%% Compute the separatrices #2
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
epsi=[0.1,0.1];
X0u2 = u(1,:)-epsi.*v1(1,:);
Tspanf2 = 1000:-0.1:0;
[~, Y2]=ode15s(@SNICeroclinic,Tspanf2,X0u2,opts, p);
%% Compute the separatrices #3
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
epsi=[1,1];
X0u2 = [-6;-2.2];
Tspanf2 = 0:0.1:1000;
[~, Y3]=ode15s(@SNICeroclinic,Tspanf2,X0u2,opts, p);
%% Compute the separatrices #4
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
epsi=[0.1,0.1];
X0u2 = u(3,:)+epsi.*v1(3,:);
Tspanf2 = 0:0.1:1000;
[~, Y4]=ode15s(@SNICeroclinic,Tspanf2,X0u2,opts, p);
%% Compute the separatrices #5
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
epsi=[0.1,0.1];
X0u2 = u(3,:)+epsi.*v2(3,:);
Tspanf2 = 1000:-0.1:0;
[~, Y5]=ode15s(@SNICeroclinic,Tspanf2,X0u2,opts, p);
%% Plot the phase plane with a saddle and a saddle node
close all
f2=figure(2);
f2.Units="centimeters";
f2.OuterPosition = [2 10 24 14];
subplot(1,2,1)
hold on; box on; grid off
set ( gca , 'FontSize' , 13 , 'fontname' , 'times');
ylabel('$y$', Interpreter='latex')
xlabel('$x$', Interpreter='latex')
% Plot separatrices
plot(Y1(1:end,1),Y1(1:end,2), 'Color','k', 'LineWidth',1.5, 'LineStyle','-')
plot(Y2(1:end,1),Y2(1:end,2), 'Color',colour.grey, 'LineWidth',1.5, 'LineStyle','-')
plot(Y3(1:end,1),Y3(1:end,2), 'Color',colour.grey, 'LineWidth',1.5, 'LineStyle','-')
plot(Y4(1:end,1),Y4(1:end,2), 'Color',colour.grey, 'LineWidth',1.5, 'LineStyle',':')
plot(Y5(1:end,1),Y5(1:end,2), 'Color',colour.grey, 'LineWidth',1.5, 'LineStyle','-')
% Plot fixed points
plot(u(1,1),u(1,2),'o','MarkerFaceColor', ...
colour.yellow,'MarkerEdgeColor','black', 'MarkerSize', 8,'LineWidth',1)
n=3;
plot(u(n,1),u(n,2),'o','MarkerFaceColor', ...
colour.pink,'MarkerEdgeColor','black', 'MarkerSize', 8,'LineWidth',1)
xlim([-6 6])
ylim([-3.5 3.5])
title('A')

%% Produce bifurcation diagram (output from AUTO 07-p)
addpath('bif_dat')
load('eq.dat')
load('eq_end.dat')
load('LC1.dat')
load('LC2.dat')
load('LC3.dat')

% Identify stability based on the analysis
for i =1:length(eq)
    Stable(i,1)=eq(i,1);
    Unstable(i,1)=eq(i,1);

    if eq(i,5) ==2
        Stable(i,2) = eq(i,3);
        Stable(i,3) = eq(i,4);
        Unstable(i,2) = NaN;
        Unstable(i,3) = NaN;

    else  
        Unstable(i,2) = eq(i,3);
        Unstable(i,3) = eq(i,4);
        Stable(i,2) = NaN;
        Stable(i,3) = NaN;

    end
end

% Location of special points (taken from AUTO 07-p)
LP=eq([325,811,1551],[1,3,4]);
HB=eq([847,1437],[1,3,4]);

subplot(1,2,2)
hold on; box on; grid off
set ( gca , 'FontSize' , 13 , 'fontname' , 'times');

ylabel('$y$', Interpreter='latex')
xlabel('$c$', Interpreter='latex')
% Plot equilibria 
plot(Stable(1:end,1),Stable(1:end,3), 'Color','k', 'LineWidth',1.5, 'LineStyle','-')
plot(Unstable(1:end,1),Unstable(1:end,3), 'Color','k', 'LineWidth',1.5, 'LineStyle',':')
plot(eq_end(1:end,1),eq_end(1:end,4), 'Color','k', 'LineWidth',1.5, 'LineStyle','-')
% Plot max and min of the oscillation amplitude
plot(LC1(1:end,1),LC1(1:end,4), 'Color','r', 'LineWidth',1.5, 'LineStyle','-')
plot(LC2(1:end,1),LC2(1:end,4), 'Color','r', 'LineWidth',1.5, 'LineStyle','-')
plot(LC3(1:end,1),LC3(1:end,4), 'Color','r', 'LineWidth',1.5, 'LineStyle','-')
plot(LC1(1:end,1),LC1(1:end,6), 'Color','b', 'LineWidth',1.5, 'LineStyle','-')
plot(LC2(1:end,1),LC2(1:end,6), 'Color','b', 'LineWidth',1.5, 'LineStyle','-')
plot(LC3(1:end,1),LC3(1:end,6), 'Color','b', 'LineWidth',1.5, 'LineStyle','-')
% Plot bifurcation points
plot(LC1(end,1),LC1(end,4),'o','MarkerFaceColor', ...
'black','MarkerEdgeColor','black', 'MarkerSize', 7,'LineWidth',1)
plot(LC1(end,1),LC1(end,6),'o','MarkerFaceColor', ...
'black','MarkerEdgeColor','black', 'MarkerSize', 7,'LineWidth',1)
plot(LC3(end,1),LC3(end,4),'o','MarkerFaceColor', ...
'black','MarkerEdgeColor','black', 'MarkerSize', 7,'LineWidth',1)
plot(LC3(end,1),LC3(end,6),'o','MarkerFaceColor', ...
'black','MarkerEdgeColor','black', 'MarkerSize', 7,'LineWidth',1)
plot(LC2(1,1),LC2(1,4),'o','MarkerFaceColor', ...
'black','MarkerEdgeColor','black', 'MarkerSize', 7,'LineWidth',1)
plot(LC2(1,1),LC2(1,6),'o','MarkerFaceColor', ...
'black','MarkerEdgeColor','black', 'MarkerSize', 7,'LineWidth',1)
plot(LC2(end,1),LC2(end,4),'o','MarkerFaceColor', ...
'white','MarkerEdgeColor','black', 'MarkerSize', 7,'LineWidth',1)
plot(LC2(end,1),LC2(end,6),'o','MarkerFaceColor', ...
'white','MarkerEdgeColor','black', 'MarkerSize', 7,'LineWidth',1)
plot(LP(:,1),LP(:,3),'^','MarkerFaceColor', ...
'black','MarkerEdgeColor','black', 'MarkerSize', 8,'LineWidth',1)
plot(HB(:,1),HB(:,3),'s','MarkerFaceColor', ...
'black','MarkerEdgeColor','black', 'MarkerSize', 8,'LineWidth',1)
xlim([-2.7 2.7])
title('B')
% saveas(f2,'SNICeroclinic.svg', 'svg')