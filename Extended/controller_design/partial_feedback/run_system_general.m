%Check the designed controller on the LCS model for the 3box example
%pathlcp needs to be in matlab path in order to run the code (we plot the
%envelope of trajectories in this example)
addpath 'C:\Users\alp1a\pathlcp\pathmexw64'

clear
clc
close all

load('controller.mat')

%extract dimension information
n = size(A,2); %dimension of state space
k = size(B,2); %dimension of input
m = size(D,2); %number of contacts

rng(0)
tend = 70;
tspan = [0 tend]; %span of a single trajectory
num_trials = 10; %number of trajectories to be tested
range = 1; %range of starting x_0 positions
tspan2 = [0 10];

for i = 1:1
y0 = range*(0.5-rand(1,n));
y0(4) = 0;
[t1,y1] = ode15s(@(t,y) sys_general(t,y,A,B,D,KK,LL,Fc,Ec), tspan, y0);
[t2,y2] = ode15s(@(t,y) sys_general(t,y,A,B,D,KK,LL,Fc,Ec), tspan2, y0);
end

%get the \lambda and u values
lambda = zeros(2, length(y1) );
lambda2 = zeros(2, length(y2) );
u = zeros(2, length(y1) );
u2 = zeros(2, length(y2) );
for i=1:length(y1)
    lambda(:,i) = pathlcp(Fc,Ec*y1(i,:)' );
    u(:,i) = KK*y1(i,:)' + LL*lambda(:,i);
end
for i=1:length(y2)
    lambda2(:,i) = pathlcp(Fc,Ec*y2(i,:)' );
    u2(:,i) = KK*y2(i,:)' + LL*lambda2(:,i);
end

LW = 2;
sz = 30;
szl = 20;
figure
subplot(2,3,1)
plot(t1,u, 'LineWidth', LW)
legend('u_1', 'u_2', 'Location','southeast','FontSize',szl)
xlabel('Time (s)')
ylabel('u_{1,2} (t)')
set(gca,'FontSize',sz);

subplot(2,3,2)
plot(t1,y1(:,1:4), 'LineWidth', LW)
legend('x_1', 'x_2', 'x_3', 'x_4', 'Location','southeast','FontSize',szl ,'NumColumns',2)
xlabel('Time (s)')
ylabel('x_{1,2,3,4} (t)')
set(gca,'FontSize',sz);

subplot(2,3,3)
plot(t1,y1(:,5:8), 'LineWidth', LW)
legend('x_5', 'x_6', 'x_7', 'x_8', 'Location','southeast','FontSize',szl ,'NumColumns',2)
xlabel('Time (s)')
ylabel('x_{5,6,7,8} (t)')
set(gca,'FontSize',sz);

subplot(2,3,4)
plot(t2,u2, 'LineWidth', LW)
legend('u_1', 'u_2', 'Location','southeast','FontSize',szl)
xlabel('Time (s)')
ylabel('u_{1,2} (t)')
set(gca,'FontSize',sz);

subplot(2,3,5)
plot(t2,y2(:,1:4), 'LineWidth', LW)
legend('x_1', 'x_2', 'x_3', 'x_4', 'Location','southeast','FontSize',szl,'NumColumns',2)
xlabel('Time (s)')
ylabel('x_{1,2,3,4} (t)')
set(gca,'FontSize',sz);

subplot(2,3,6)
plot(t2,y2(:,5:8), 'LineWidth', LW)
legend('x_5', 'x_6', 'x_7', 'x_8', 'Location','southeast','FontSize',szl ,'NumColumns',2)
xlabel('Time (s)')
ylabel('x_{5,6,7,8} (t)')
set(gca,'FontSize',sz);

% rng shuffle
% subplot(4,3,[7 8 9])
% for i = 1:num_trials
%     i
%     y0 = range*(0.5-rand(1,n));
%     [time{i},y{i}] = ode45(@(t,y) sys_general(t,y,A,B,D,KK,LL,Fc,Ec), tspan, y0);
%     plot(time{i},y{i}, 'LineWidth', LW, 'Color', [0.5, 0.5, 0.5] )
%     hold on
% end
% 
% xlabel('Time (s)')
% ylabel('\{x(t)\}_i')

% %plot the envelope
% prec = 0.01;
% t = 0:prec:tend;
% y_max = zeros(1,length(t)-1);
% y_min = zeros(1,length(t)-1);
% for i = 1:length(t)-1  %go through bins 
%     %go through each trajectory
%     for j = 1:num_trials
%         for k = 1:length(time{j})
%             if time{j}(k) >= t(i) && time{j}(k) <= t(i+1)
%                 if max ( y{j}(k,:) ) >= y_max(i)
%                    y_max(i) = max ( y{j}(k,:) ); 
%                 end
%                 if min ( y{j}(k,:) ) <= y_min(i)
%                     y_min(i) = min ( y{j}(k,:) );
%                 end
%             end
%         end
%     end
% end
% 
% subplot(4,3,[7 8 9])
% plot(0:prec:tend-prec,y_max,'LineWidth',1,'Color',[0.5,0.5,0.5])
% hold on
% plot(0:prec:tend-prec,y_min,'LineWidth',1,'Color',[0.5,0.5,0.5])
% hold on
% area(0:prec:tend-prec,y_max,'FaceColor',[0.5,0.5,0.5],'EdgeColor', [0.5,0.5,0.5],'FaceAlpha',0.3);
% hold on
% area(0:prec:tend-prec,y_min,'FaceColor',[0.5,0.5,0.5],'EdgeColor', [0.5,0.5,0.5],'FaceAlpha',0.3);
% 
% xlabel('Time (s)')
% ylabel('\{x(t)\}_i')