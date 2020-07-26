clear all
clc
addpath 'C:\Users\alp1a\pathlcp\pathmexw64'

tactile_feedback = 0; %1 to use tactile feedback, 0 for LQR

%parameters
g = 9.81;
l1 = 0.5;
l2 = 1;
m1 = 0.5;
m2 = 1;
a = m1*(l1^2) + m2*(l1^2);
b = m2*(l2^2);
c = m2*l1*l2;
d = g*m1*l1 + g*m2*l1;
e = g*m2*l2;
angle_cons = 0.2;
d1 = angle_cons;
d2 = -angle_cons;

Ma = [a + b + 2*c b+c; b+c b];
Jd = [-1 1; 0 0];

%problem data
%state matrix, input matrix, contact matrix
A =  [zeros(2) eye(2); g/l1 -(g*m2)/(l1*m1) 0 0; -g/l1 g*(l1*m1+l1*m2+l2*m2)/(l1*l2*m1) 0 0]; 
B = [zeros(2,1); -(l1+l2)/(l2*m1*(l1^2)); ( (m1* (l1^2 ) )+m2*(l1+l2)^2 )/((l1^2)*(l2^2)*m1*m2)];
D = [zeros(2,2); inv(Ma)*Jd]; Ec = [-1 0 0 0; 1 0 0 0]; Fc = 1*eye(2); w = [d1; -d2];
 
%LQR
KK = -lqr(A,B,100*eye(4),1); 
LL = [ 0  0];
%discrete euler
%rng(1) %example where nonlinear works and linear fails
rng(1)
x_save = [];
lamh = [];
%discrete euler (nonlinear)
    dt = 0.001; %stepsize
    t = 2.2; %simulation time
    x = zeros(4,round(t/dt) + 2);
    u = zeros(1,round(t/dt) + 2);
    x(:,1) = [0.0416428766398130;-0.0609530023780252;1.61508671851095;-2.65346218839446];
    for i = 1:(t/dt)+1
        x1 = x(1,i); x2 = x(2,i); x3 = x(3,i); x4 = x(4,i);
        phi1 = d1 - x1;
        phi2 = x1 - d2;
        lambda = pathlcp(Fc,[phi1;phi2]);
        lamh = [lamh lambda];
        %lambda = [0;0];
        if max(lambda) >= 10^6
            flag2 = 1;
            counter2 = counter2 + 1;
            x(:)=10;
            break
        end
        M = [a+b+2*c*cos(x2) b+c*cos(x2);b+c*cos(x2) b];
        C = [-c*sin(x2)*x4 -c*sin(x2)*(x3+x4); c*sin(x2)*x3 0];
        G = [-d*sin(x1) - e*sin(x1+x2); e*sin(x1+x2)]; 
        Bb =[0;1];
        J = [-1 1; 0 0];
        Minv = inv(M); B = [zeros(2,1); Minv*Bb];
        x(3:4,i+1) = x(3:4,i) + (-Minv*G - Minv*C*[x3;x4] + Minv*J*lambda)*dt;
        x(1,i+1) = x(1,i) + x(3,i)*dt;
        x(2,i+1) = x(2,i) + x(4,i)*dt;
        %controllers
        cont = [x1; x2; x3; x4];
        x(:,i+1) = x(:,i+1) + (B*KK)*cont*dt + (B*LL)*lambda*dt; 
        u(i) = (KK)*cont + (LL)*lambda;
    end
    
LW = 4;
sz = 30;
szl = 20;
figure
subplot(2,2,1)
plot([0:dt:t],x(:,1:end-1),'LineWidth',LW)
legend('x1','x2','x3','x4', 'Location','southwest','FontSize',szl,'NumColumns',2)
xlabel('Time (s)')
ylabel('x(t)')
title('LQR','FontSize',sz)
set(gca,'FontSize',sz);

subplot(2,2,2)
plot([0:dt:t],u(:,1:end-1),'LineWidth',LW)
legend('u','FontSize',szl)
xlabel('Time (s)')
ylabel('u(t)')
title('LQR','FontSize',sz)
set(gca,'FontSize',sz);

%Contact-Aware
KK = [726.1180  438.3622  283.9890  174.2047];
LL = [0.3859   -0.3852];
%discrete euler
x_save = [];
lamh = [];
%discrete euler (nonlinear)
    dt = 0.001; %stepsize
    t = 5; %simulation time
    x = zeros(4,round(t/dt) + 2);
    u = zeros(1,round(t/dt) + 2);
    x(:,1) = [0.0416428766398130;-0.0609530023780252;1.61508671851095;-2.65346218839446];
    for i = 1:(t/dt)+1
        x1 = x(1,i); x2 = x(2,i); x3 = x(3,i); x4 = x(4,i);
        phi1 = d1 - x1;
        phi2 = x1 - d2;
        lambda = pathlcp(Fc,[phi1;phi2]);
        lamh = [lamh lambda];
        %lambda = [0;0];
        if max(lambda) >= 10^6
            flag2 = 1;
            counter2 = counter2 + 1;
            x(:)=10;
            break
        end
        M = [a+b+2*c*cos(x2) b+c*cos(x2);b+c*cos(x2) b];
        C = [-c*sin(x2)*x4 -c*sin(x2)*(x3+x4); c*sin(x2)*x3 0];
        G = [-d*sin(x1) - e*sin(x1+x2); e*sin(x1+x2)]; 
        Bb =[0;1];
        J = [-1 1; 0 0];
        Minv = inv(M); B = [zeros(2,1); Minv*Bb];
        x(3:4,i+1) = x(3:4,i) + (-Minv*G - Minv*C*[x3;x4] + Minv*J*lambda)*dt;
        x(1,i+1) = x(1,i) + x(3,i)*dt;
        x(2,i+1) = x(2,i) + x(4,i)*dt;
        %controllers
        cont = [x1; x2; x3; x4];
        x(:,i+1) = x(:,i+1) + (B*KK)*cont*dt + (B*LL)*lambda*dt; 
        u(i) = (KK)*cont + (LL)*lambda;
    end
    
subplot(2,2,3)
plot([0:dt:t],x(:,1:end-1),'LineWidth',LW)
legend('x1','x2','x3','x4', 'Location','southeast','FontSize',szl,'NumColumns',2)
xlabel('Time (s)')
ylabel('x(t)')
title('Contact-Aware Policy','FontSize',sz)
set(gca,'FontSize',sz);
xlim([0 5])

subplot(2,2,4)
plot([0:dt:t],u(:,1:end-1),'LineWidth',LW)
legend('u','FontSize',szl, 'Location','southeast')
xlabel('Time (s)')
ylabel('u(t)')
title('Contact-Aware Policy','FontSize',sz)
set(gca,'FontSize',sz);
xlim([0 5])

%plot all trajectories for tactile

% %discrete euler
% rng(0)
% subplot(3,2,[5 6])
% number_of_trials = 100; %number of initial conditions
% for k = 1:number_of_trials
%     k %to see the progress
%     flag2 = 0; %is 1 if \lambda exceeds a certain value
%     dt = 0.001; %stepsize
%     t = 5; %simulation time
%     x = zeros(4,t/dt + 2); %x values
%     x(1,1) = 0*(rand(1)-0.5); x(2,1) = 0*(rand(1)-0.5); x(3,1) = 0.1*(rand(1)-0.5); x(4,1) = 0.1*(rand(1)-0.5);
%     %Run the simulaton
%     for i = 1:(t/dt)+1
%         x1 = x(1,i); x2 = x(2,i); x3 = x(3,i); x4 = x(4,i);
%         phi1 = d1 - x1;
%         phi2 = x1 - d2;
%         lambda = pathlcp(Fc,[phi1;phi2]);
%         lamh = [lamh lambda];
%         %lambda = [0;0];
%         if max(lambda) >= 10^6
%             flag2 = 1;
%             counter2 = counter2 + 1;
%             x(:)=10;
%             break
%         end
%         M = [a+b+2*c*cos(x2) b+c*cos(x2);b+c*cos(x2) b];
%         C = [-c*sin(x2)*x4 -c*sin(x2)*(x3+x4); c*sin(x2)*x3 0];
%         G = [-d*sin(x1) - e*sin(x1+x2); e*sin(x1+x2)]; 
%         Bb =[0;1];
%         J = [-1 1; 0 0];
%         Minv = inv(M); B = [zeros(2,1); Minv*Bb];
%         x(3:4,i+1) = x(3:4,i) + (-Minv*G - Minv*C*[x3;x4] + Minv*J*lambda)*dt;
%         x(1,i+1) = x(1,i) + x(3,i)*dt;
%         x(2,i+1) = x(2,i) + x(4,i)*dt;
%         %controllers
%         cont = [x1; x2; x3; x4];
%         x(:,i+1) = x(:,i+1) + (B*KK)*cont*dt + (B*LL)*lambda*dt; 
%         u(i) = (KK)*cont + (LL)*lambda;
%     end
%     if abs(max(x(:,i))) < 0.1
%         plot([0:dt:t],x(:,1:end-1),'LineWidth',LW, 'Color', [0.5, 0.5, 0.5])
%         hold on
%     end
% end
% 
% xlabel('Time (s)')
% ylabel('\{x(t)\}_i')