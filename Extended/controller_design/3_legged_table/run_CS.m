clear all
clc
close all

load('controller.mat')
%addpath PATHLCP

%extract dimension information
n = size(A,2); %dimension of state space
k = size(B,2); %dimension of input
m = size(D,2); %number of contacts

tspan = [0 2]; %span of a single trajectory

LW = 4;
y0(1) = -1;
y0(2) = 20;
[t,y] = ode15s(@(t,y) sys_affine(t,y,A,B,D,KK,LL,m,Fc,Ec,c,kappa,H,k), tspan, y0);

sz = 30;
figure
subplot(1,2,1)
plot(t,y(:,2), 'LineWidth', LW)
%legend('\tau', 'Location','southeast','FontSize',20)
xlabel('Time (s)')
ylabel('\tau(t)')
set(gca,'FontSize',sz);

subplot(1,2,2)
plot(t,y(:,1), 'LineWidth', LW)
%legend('x', 'Location','northeast','FontSize',20 ,'NumColumns',2)
xlabel('Time (s)')
ylabel('x(t)')
set(gca,'FontSize',sz);