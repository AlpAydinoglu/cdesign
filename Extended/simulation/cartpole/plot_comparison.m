clear all
clc
addpath 'C:\Users\alp1a\pathlcp\pathmexw64'

%parameters
g = 9.81;
mp = 0.1;
mc = 1;
len = 0.5;
alpha = pi/2;
alpha2 = pi/2;
r1 = 0;  
r2 = 0;
d1 = 0.1;
d2 = -0.1;
w1 = d1*sin(alpha2) + len*cos(alpha2) - r1*cos(alpha2);
w2 = -d2*sin(alpha) + len*cos(alpha) - r2*cos(alpha);

%problem data (FOR 2 WALLS)
%state matrix, input matrix, contact matrix
A =  [zeros(2) eye(2); 0 g*mp/mc 0 0; 0 g*(mc+mp)/(len*mc) 0 0]; %state matrix
B = [0; 0; 1/mc; 1/(len*mc)]; %input matrix
D = [zeros(3,2); sin(alpha2)/(len*mp) -sin(alpha)/(len*mp)]; %contact matrix
%complementarity constraints
Ec = [-1*sin(alpha2) len*sin(alpha2) 0 0; 1*sin(alpha) -len*sin(alpha) 0 0]; 
Fc = 0.1*eye(2);  %STIFFNESS PARAMETER HERE
w = [w1;w2];
 
%LQR
KK = -lqr(A,B,100*eye(4),1); 
LL = [ 0   0];
rng(0)
%discrete euler
number_of_trials = 1; %number of initial conditions
counter = 0; %counts number of times the system is unstable
counter2 = 0;
    flag = 0; %1 if contact force exceeds a certain value
    dt = 0.01; %stepsize
    t = 2.35; %simulation time
    x = zeros(4,round(t/dt) + 2); %holds state values
    u = zeros(1,round(t/dt) + 2);
    %initial conditions
    x(1,1) = 0.2*(rand(1)-0.5); x(2,1) = 0+pi; x(3,1) = 2*(rand(1)-0.5); x(4,1) = 2*(rand(1)-0.5);
    %Run the simulaton
    for i = 1:(t/dt)+1
        %x1,x2,x3,x4 (\bar{x2} = x2 - pi)
        x1 = x(1,i); x2 = x(2,i)-pi; x3 = x(3,i); x4 = x(4,i);
        
        %calculate coordinates of the pole
        p = [x1 - len*sin(x2); len*cos(x2)];
        
        %Gap functions phi1, phi2
        phi1 = d1*sin(alpha) - p(1)*sin(alpha) + p(2)*cos(alpha)-r1*cos(alpha);
        phi2 = p(1)*sin(alpha) - d2*sin(alpha) + p(2)*cos(alpha) - r2*cos(alpha);
        
        %calculate the contact force lambda
        lambda = pathlcp(Fc,[phi1;phi2]);
        
        %Calculate M,C,G,J,B
        x2 = x2 + pi;
        M = [mc+mp mp*len*cos(x2); mp*len*cos(x2) mp*(len^2)];
        C = [0 -mp*len*x4*sin(x2); 0 0];
        G = [0; mp*g*len*sin(x2)]; Bb =[1;0];
        J = [-sin(alpha) sin(alpha); len*sin(alpha)-len*cos(alpha)*x2 -len*sin(alpha)-len*cos(alpha)*x2];
        Minv = inv(M); B = [zeros(2,1); Minv*Bb];
        
        %Iterate the autonomous part of the dynamics
        x(3:4,i+1) = x(3:4,i) + (-Minv*G - Minv*C*[x3;x4] + Minv*J*lambda)*dt;
        x(1,i+1) = x(1,i) + x(3,i)*dt;
        x(2,i+1) = x(2,i) + x(4,i)*dt;
        
        %Add the control action
        cont = [x1; x2-pi; x3; x4];
        x(:,i+1) = x(:,i+1) + (B*KK)*cont*dt + (B*LL)*lambda*dt; 
        u(i) = (KK)*cont + (LL)*lambda;
    end
    
x(2,:) = x(2,:)-pi;

sz = 30;
szl = 20;
LW = 4;
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
xlabel('Time (s)')
ylabel('u(t)')
title('LQR','FontSize',sz)
set(gca,'FontSize',sz);

%tactile feedback
KK = [3.6931  -46.7030    3.3885   -5.7147];
LL = [-13.9797   13.9767];
rng(0)

%discrete euler
number_of_trials = 1; %number of initial conditions
counter = 0; %counts number of times the system is unstable
counter2 = 0;
    flag = 0; %1 if contact force exceeds a certain value
    dt = 0.01; %stepsize
    t = 6; %simulation time
    x = zeros(4,round(t/dt) + 2); %holds state values
    u = zeros(1,round(t/dt) + 2);
    %initial conditions
    x(1,1) = 0.2*(rand(1)-0.5); x(2,1) = 0+pi; x(3,1) = 2*(rand(1)-0.5); x(4,1) = 2*(rand(1)-0.5);
    %Run the simulaton
    for i = 1:(t/dt)+1
        %x1,x2,x3,x4 (\bar{x2} = x2 - pi)
        x1 = x(1,i); x2 = x(2,i)-pi; x3 = x(3,i); x4 = x(4,i);
        
        %calculate coordinates of the pole
        p = [x1 - len*sin(x2); len*cos(x2)];
        
        %Gap functions phi1, phi2
        phi1 = d1*sin(alpha) - p(1)*sin(alpha) + p(2)*cos(alpha)-r1*cos(alpha);
        phi2 = p(1)*sin(alpha) - d2*sin(alpha) + p(2)*cos(alpha) - r2*cos(alpha);
        
        %calculate the contact force lambda
        lambda = pathlcp(Fc,[phi1;phi2]);
        
        %Calculate M,C,G,J,B
        x2 = x2 + pi;
        M = [mc+mp mp*len*cos(x2); mp*len*cos(x2) mp*(len^2)];
        C = [0 -mp*len*x4*sin(x2); 0 0];
        G = [0; mp*g*len*sin(x2)]; Bb =[1;0];
        J = [-sin(alpha) sin(alpha); len*sin(alpha)-len*cos(alpha)*x2 -len*sin(alpha)-len*cos(alpha)*x2];
        Minv = inv(M); B = [zeros(2,1); Minv*Bb];
        
        %Iterate the autonomous part of the dynamics
        x(3:4,i+1) = x(3:4,i) + (-Minv*G - Minv*C*[x3;x4] + Minv*J*lambda)*dt;
        x(1,i+1) = x(1,i) + x(3,i)*dt;
        x(2,i+1) = x(2,i) + x(4,i)*dt;
        
        %Add the control action
        cont = [x1; x2-pi; x3; x4];
        x(:,i+1) = x(:,i+1) + (B*KK)*cont*dt + (B*LL)*lambda*dt; 
        u(i) = (KK)*cont + (LL)*lambda;
    end

x(2,:) = x(2,:)-pi;
    
subplot(2,2,3)
plot([0:dt:t],x(:,1:end-1),'LineWidth',LW)
legend('x1','x2','x3','x4', 'Location','southeast','FontSize',szl,'NumColumns',2)
xlabel('Time (s)')
ylabel('x(t)')
title('Contact-Aware Policy','FontSize',sz)
set(gca,'FontSize',sz);

subplot(2,2,4)
plot([0:dt:t],u(:,1:end-1),'LineWidth',LW)
xlabel('Time (s)')
ylabel('u(t)','FontSize',szl)
title('Contact-Aware Policy','FontSize',sz)
set(gca,'FontSize',sz);


%plot all trajectories for tactile

%discrete euler
rng(0)
% subplot(3,2,[5 6])
number_of_trials = 1; %number of initial conditions
counter = 0; %counts number of times the system is unstable
counter2 = 0;
for k = 1:number_of_trials
    k
    flag = 0; %1 if contact force exceeds a certain value
    dt = 0.01; %stepsize
    t = 6; %simulation time
    x = zeros(4,round(t/dt) + 2); %holds state values
    u = zeros(1,round(t/dt) + 2);
    %initial conditions
    x(1,1) = 0.2*(rand(1)-0.5); x(2,1) = 0+pi; x(3,1) = 2*(rand(1)-0.5); x(4,1) = 2*(rand(1)-0.5);
    %Run the simulaton
    for i = 1:(t/dt)+1
        %x1,x2,x3,x4 (\bar{x2} = x2 - pi)
        x1 = x(1,i); x2 = x(2,i)-pi; x3 = x(3,i); x4 = x(4,i);
        
        %calculate coordinates of the pole
        p = [x1 - len*sin(x2); len*cos(x2)];
        
        %Gap functions phi1, phi2
        phi1 = d1*sin(alpha) - p(1)*sin(alpha) + p(2)*cos(alpha)-r1*cos(alpha);
        phi2 = p(1)*sin(alpha) - d2*sin(alpha) + p(2)*cos(alpha) - r2*cos(alpha);
        
        %calculate the contact force lambda
        lambda = pathlcp(Fc,[phi1;phi2]);
        
        %Calculate M,C,G,J,B
        x2 = x2 + pi;
        M = [mc+mp mp*len*cos(x2); mp*len*cos(x2) mp*(len^2)];
        C = [0 -mp*len*x4*sin(x2); 0 0];
        G = [0; mp*g*len*sin(x2)]; Bb =[1;0];
        J = [-sin(alpha) sin(alpha); len*sin(alpha)-len*cos(alpha)*x2 -len*sin(alpha)-len*cos(alpha)*x2];
        Minv = inv(M); B = [zeros(2,1); Minv*Bb];
        
        %Iterate the autonomous part of the dynamics
        x(3:4,i+1) = x(3:4,i) + (-Minv*G - Minv*C*[x3;x4] + Minv*J*lambda)*dt;
        x(1,i+1) = x(1,i) + x(3,i)*dt;
        x(2,i+1) = x(2,i) + x(4,i)*dt;
        
        %Add the control action
        cont = [x1; x2-pi; x3; x4];
        x(:,i+1) = x(:,i+1) + (B*KK)*cont*dt + (B*LL)*lambda*dt; 
        u(i) = (KK)*cont + (LL)*lambda;
    end
    x(2,:) = x(2,:)-pi;
%    plot([0:dt:t],x(:,1:end-1),'LineWidth',LW, 'Color', [0.5, 0.5, 0.5])
%    hold on
end

% xlabel('Time (s)')
% ylabel('\{x(t)\}_i')