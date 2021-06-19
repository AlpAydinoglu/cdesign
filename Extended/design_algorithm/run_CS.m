clear all
clc
close all

%use this code to evaluate the system only if the x-trajectory is unique.
%if x-trajectory is not unique, you are dealing with a differential
%inclusion be cautious, this would approximate a single solution (there
%could be jumps in the solution so ODE45 might not be the best tool)

load('controller.mat') %gets controller and system matrices
%addpath
%here pathlcp should be in the path

%extract dimension information
n = size(A,2); %dimension of state space
k = size(B,2); %dimension of input
m = size(D,2); %number of contacts

tspan = []; %span of a single trajectory (tspan has the form [ a b ] where t_0 = a and t_f = b )
y0 = []; %initial condition
[t,y] = ode15s(@(t,y) sys_affine(t,y,A,B,D,KK,LL,m,Fc,Ec,c,k), tspan, y0); %evaluate the system