clear
clc
close all

load('controller.mat')
addpath 'C:\Users\Alp\Desktop\linear complementarity systems\codes\path'

%extract dimension information
n = size(A,2); %dimension of state space
k = size(B,2); %dimension of input
m = size(D,2); %number of contacts

tend = 20;
tspan = [0 tend]; %span of a single trajectory
num_trials = 2; %number of trajectories to be tested
%KK = zeros(1,2); LL = zeros(1,3); %for trying the case where there is no
range = 0.2; %range of starting x_0 positions
%feedback
%simulate and plot
%fornow
G{1} = zeros(n); G{2} = zeros(n); G{3} = zeros(n); 
JJ{1} = zeros(k,n); JJ{2} = zeros(k,n); JJ{3} = zeros(k,n);

for i = 1:num_trials
y0 = range*(0.5-rand(1,n));
[time{i},y{i}] = ode45(@(t,y) sys_affine(t,y,A,B,D,KK,LL,G,JJ,m,Fc,Ec,w), tspan, y0);
end

%plot the envelope
prec = 0.1;
t = 0:prec:tend;
y_max = zeros(1,length(t)-1);
y_min = zeros(1,length(t)-1);
for i = 1:length(t)-1  %go through bins 
    %go through each trajectory
    for j = 1:num_trials
        for k = 1:length(time{j})
            if time{j}(k) >= t(i) && time{j}(k) <= t(i+1)
                if max ( y{j}(k,:) ) >= y_max(i)
                   y_max(i) = max ( y{j}(k,:) ); 
                end
                if min ( y{j}(k,:) ) <= y_min(i)
                    y_min(i) = min ( y{j}(k,:) );
                end
            end
        end
    end
end

figure
plot(0:0.1:tend-prec,y_max,'LineWidth',1,'Color',[0.5,0.5,0.5])
hold on
plot(0:0.1:tend-prec,y_min,'LineWidth',1,'Color',[0.5,0.5,0.5])
hold on
area(0:0.1:tend-prec,y_max,'FaceColor',[0.5,0.5,0.5],'EdgeColor', [0.5,0.5,0.5],'FaceAlpha',0.3);
hold on
area(0:0.1:tend-prec,y_min,'FaceColor',[0.5,0.5,0.5],'EdgeColor', [0.5,0.5,0.5],'FaceAlpha',0.3);

for i = 1:3
hold on
plot(time{i},y{i},'Color',[0.5,0.5,0.5],'LineWidth',0.1) 
end

xlabel('Time (s)')
ylabel('\{x(t)\}_i')