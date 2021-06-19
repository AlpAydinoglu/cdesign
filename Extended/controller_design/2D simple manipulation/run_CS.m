clear all
clc
close all

load('controller.mat')
%addpath pathlcp

%extract dimension information
n = size(A,2); %dimension of state space
k = size(B,2); %dimension of input
m = size(D,2); %number of contacts

tspan = [0 5]; %span of a single trajectory

y0(1) = 5; y0(2) = 0; y0(3) = 10; 
tau_zero = KK*[y0(1);y0(2);y0(3)] + LL*pathlcp(Fc,Ec*[y0(1);y0(2);y0(3)]+c);
y0(4) = tau_zero(1); y0(5) = tau_zero(2);
[t,y] = ode15s(@(t,y) sys_affine(t,y,A,B,D,KK,LL,m,Fc,Ec,c,kappa,H,k), tspan, y0);
LW = 4;
sz = 30;
figure
subplot(1,2,1)
plot(t,y(:,4:5), 'LineWidth', LW)
legend('\tau_{1}', '\tau_{2}', 'Location','southeast','FontSize',30)
xlabel('Time (s)')
ylabel('\tau(t)')
set(gca,'FontSize',sz);
xlim([0 5])

subplot(1,2,2)
plot(t,y(:,1:3), 'LineWidth', LW)
legend('x_1','x_2','x_3', 'Location','northeast','FontSize',30)
xlabel('Time (s)')
ylabel('x(t)')
set(gca,'FontSize',sz);
xlim([0 5])

% subplot(2,2,[3 4])
% num_trials = 10;
% range = 10; %range of starting x_0 positions
% %simulate and plot
% for i = 1:num_trials
% y0 = range*(0.5-rand(1,n+k));
% [t,y] = ode45(@(t,y) sys_affine(t,y,A,B,D,KK,LL,m,Fc,Ec,c,kappa,H,k), tspan, y0);
% plot(t,y(:,1:3),'LineWidth',LW, 'Color', [0.5, 0.5, 0.5])
% hold on
% end
% xlabel('Time (s)')
% ylabel('\{x(t)\}_i')
% 
% % plot(t,y(:,1:3),'LineWidth',2)
% % legend('x', 'x_l', 'x_r')
% % xlabel('Time (s)')
% % %ylabel('\{x_1(t), x_2(t) \}_i')
% % figure
% % plot(t,y(:,4:5),'LineWidth',2)
% 
% 
lam_val = [];
for i = 1:length(y)
    x_val = y(i,1:end-k);
    tau_val = y(i,end-k+1:end);
    lam_val = [lam_val pathlcp(Fc,Ec*x_val'+c+H*tau_val')];
end

figure
plot(t, lam_val(3:5,:), 'LineWidth',1)
legend('\gamma', '\lambda_{t}^+', '\lambda_{t}^-')