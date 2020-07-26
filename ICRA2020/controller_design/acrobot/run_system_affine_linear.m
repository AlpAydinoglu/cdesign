%Check the designed controller on the LCS model for the acrobot example
%pathlcp needs to be in matlab path in order to run the code

clear
clc
close all

%extract system and controller parameters
load('controller.mat')

%extract dimension information
n = size(A,2); %dimension of state space
k = size(B,2); %dimension of input
m = size(D,2); %number of contacts

tspan = [0 10]; %span of a single trajectory
num_trials = 1; %number of trajectories to be tested
range = 100; %range of starting x_0 positions
%feedback
%simulate and plot

for i = 1:num_trials
y0 = range*(0.5-rand(1,n));
y0(3)=0;
y0(4)=0;
[t,y] = ode45(@(t,y) sys_affine(t,y,A,B,D,KK,LL,Fc,Ec,w), tspan, y0);
hold on
plot(t,y,'Color',[0.5,0.5,0.5],'LineWidth',0.5)
end
xlabel('Time (s)')
ylabel('\{x(t)\}_i')