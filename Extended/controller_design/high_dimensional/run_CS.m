clear all
clc
close all

load('controller.mat')
addpath 'C:\Users\alp1a\pathlcp\pathmexw64'

%extract dimension information
n = size(A,2); %dimension of state space
k = size(B,2); %dimension of input
m = size(D,2); %number of contacts

tspan = [0 3]; %span of a single trajectory
y0 = [ -0.0932   -0.2165    0.1105    0.1958    0.0242    0.1078    0.0128   -0.1052];
[t,y] = ode45(@(t,y) sys_affine(t,y,A,B,D,KK,LL,m,Fc,Ec,c,k), tspan, y0);
u = zeros(6, length(y) );
lambda = zeros(10, length(y) );
for i=1:length(y)
    lambda(:,i) = pathlcp(Fc,Ec*y(i,:)' );
    u(:,i) = KK*y(i,:)' + LL*lambda(:,i);
end

LW = 4;
sz = 30;
szl = 20;
lmt = 500;
figure
subplot(1,2,1)
plot(t(1:lmt),u(:,1:lmt), 'LineWidth', LW)
legend('u_1', 'u_2', 'u_3', 'u_4', 'u_5', 'u_6', 'Location','southeast','FontSize',szl)
xlabel('Time (s)')
ylabel('u(t)')
set(gca,'FontSize',sz);

subplot(1,2,2)
plot(t(1:lmt),y(1:lmt,:), 'LineWidth', LW)
legend('x_1','x_2','x_3','x_4','x_5','x_6','x_7','x_8', 'Location','southeast','FontSize',szl,'NumColumns',2)
xlabel('Time (s)')
ylabel('x(t)')
set(gca,'FontSize',sz);


% num_trials = 5; %number of trajectories to be tested
% range = 1; %range of starting x_0 positions
% %simulate and plot
% rng shuffle
% subplot(2,2,[3 4])
% for i = 1:num_trials
% y0 = range*(0.5-rand(1,n));
% [t,y] = ode45(@(t,y) sys_affine(t,y,A,B,D,KK,LL,m,Fc,Ec,c,k), [0 5], y0);
% hold on
% plot(t,y,'LineWidth',LW, 'Color', [0.5, 0.5, 0.5])
% end
% xlabel('Time (s)')
% ylabel('\{x(t)\}_i')