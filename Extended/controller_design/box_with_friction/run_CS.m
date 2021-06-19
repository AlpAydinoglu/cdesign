clear all
clc
close all

load('controller.mat')
%addpath pathlcp

%extract dimension information
n = size(A,2); %dimension of state space
k = size(B,2); %dimension of input
m = size(D,2); %number of contacts

tspan = [0 4]; %span of a single trajectory

LW = 4;
y0(1) = 1;
y0(2) = -15;
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

% subplot(2,2,[3 4])
% num_trials = 10; %number of trajectories to be tested
% range = 50; %range of starting x_0 positions
% %simulate and plot
% for i = 1:num_trials
% y0 = range*(0.5-rand(1,n+k));
% [t,y] = ode45(@(t,y) sys_affine(t,y,A,B,D,KK,LL,m,Fc,Ec,c,kappa,H,k), tspan, y0);
% plot(t,y(:,1),'LineWidth',LW, 'Color', [0.5, 0.5, 0.5])
% hold on
% end
% xlabel('Time (s)')
% ylabel('\{x(t)\}_i')
% 
lam_val = [];
for i = 1:length(y)
    x_val = y(i,1);
    tau_val = y(i,2);
    lam_val = [lam_val pathlcp(Fc,Ec*x_val+c+H*tau_val)];
end

figure
plot(t, lam_val, 'LineWidth',1)
legend('\gamma', '\lambda_{t}^+', '\lambda_{t}^-')
% % figure
% % plot(t, lam_val(2,:)-lam_val(3,:), 'LineWidth',1)
% 
% %mu = 0.1
% P1 = 229.3534;
% P2 = 0.54183;
% P3 = -0.025006;
% P4 = 3.9809;
% P5 = -0.0089216;
% P6 = 0.10925;
% W = [0 1 -1];
% 
% %mu = 5 
% % W = [0 1 -1];
% % P1 = 219.2651;
% % P2 = 0.67718;
% % P3 = -0.03446;
% % P4 = 4.0866;
% % P5 = -0.012859;
% % P6 = 0.1175;
% 
% %mu = 5 \kappa = 1
% % W = [0 1 -1];
% % P1 = 236.4161;
% % P2 = 3.1351;
% % P3 = 6.5684;
% % P4 = 5.8212;
% % P5 = 6.711;
% % P6 = 6.975;
% 
% lam_val = [];
% V = [];
% V_dot = [];
% for i = 1:length(y)
%     x_val = y(i,1);
%     tau_val = y(i,2);
%     x_basis = x_val;
%     tau_basis = tau_val;
%     lam_basis = pathlcp(Fc,Ec*x_val+c+H*tau_val);
%     V_s = x_basis' * P1 * x_basis + 2 * x_basis' * P2 * W * lam_basis ...
%     + lam_basis' * W' * P3 * W * lam_basis + 2 * x_basis' * P4 * tau_basis ...
%     + 2 * lam_basis' * W' * P5 * tau_basis + tau_basis' * P6 * tau_basis;
%     V = [V V_s]; 
% end
% 
% % figure
% % plot(t, V, 'LineWidth',1)
% % legend('Lyapunov function')
% % 
% % Vdot = [diff(V) 0];
% % figure
% % plot(t, Vdot, 'LineWidth',1)
% % legend('Lyapunov function derivative')