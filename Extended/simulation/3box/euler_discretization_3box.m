clear all
clc

%parameters of the system
mc = 1;
m1 = 1;
m2 = 1;
m3 = 1;
g = 9.81;
mp = 1.5;
len = 0.5;

%problem data
%state matrix, input matrix, contact matrix
A =  [zeros(4) eye(4); zeros(4) zeros(4)]; A(8,4) = g*(m1+mp)/(len*m1); A(5,4) = g*mp/m1;
B = [zeros(4,2); 1/m1 0; 0 0; 0 1/m3; 1/(m1*len) 0];
D = [zeros(4,2); -1/m1 0; 1/m2 -1/m2; 0 1/m3; -1/(len*m1) 0];
%linear complementarity matrices
Ec = [-1 1 0 0 0 0 0 0; 0 -1 1 0 0 0 0 0]; Fc = 0.01*eye(2);

%extract dimension information
n = size(A,2); %dimension of state space
k = size(B,2); %dimension of input
m = size(D,2); %number of contacts

%gain matrices (precalculated)
KK = [5.2182    0.0000    0.6848 -353.2109    8.9340    0.0000   -3.1900  -55.1658; 6.8662   -0.0000   -6.5382  -51.2013    2.3006   -0.0000   -2.7591  -10.0627];
LL = [-4.0164   -0.2520; -0.2520    6.4276];

rng(1) %set rng for consistency
x_save = []; %save x values
lamh = []; %save \lam values
%discrete euler
number_of_trials = 1; %number of total trials to do
counter = 0; %counts number of times the system is unstable
counter2 = 0; %internal counter
for k = 1:number_of_trials
    k %shows the current trial
    flag = 0; %1 if contact force exceeds a certain value
    dt = 0.01; %stepsize
    t = 50; %simulation time
    x = zeros(4,t/dt + 2); %holds state values
    q_val = zeros(4,t/dt + 2);
    %initial conditions
    x(1:8,1) = 0.1*rand(8,1); x(4,1) = -0.5; x(1,1) = 0.1; x(2,1) = -0.1; x(3,1) = -0.1;
    q_val(1:8,1) = x(1:8,1); 
    q_val(2,1) = x(4,1); q_val(3,1) = x(2,1); q_val(4,1) = x(3,1);
    x(4,1) = x(4,1) + pi;
    q_val(2,1) = q_val(2,1) + pi;
    %Run the simulaton
    for i = 1:(t/dt)+1
        %change basis q1 = x1, q2 = \theta, q3 = x2, q4 = x3
        q1 = x(1,i); q2 = x(4,i); q3 = x(2,i); q4 = x(3,i);
        q5 = x(5,i); q6 = x(8,i); q7 = x(6,i); q8 = x(7,i);
        q2 = q2 - pi;
            
        %Gap functions phi1, phi2
        phi1 = q3 - q1;
        phi2 = q4 - q3;
        
        %calculate the contact force lambda
        lambda = pathlcp(Fc,[phi1;phi2]);
        %lamh = [lamh lambda]; %to hold values \lam if preferred
        
        %break if contact force exceeds 10^7
        if max(lambda) >= 10^7
            flag = 1;
            counter = counter + 1;
            x(:) = 10;
            break
        end
        
        %Calculate M,C,G,J,B,u
        q2 = q2 + pi;
        M = [mc+mp mp*len*cos(q2); mp*len*cos(q2) mp*(len^2)];
        C = [0 -mp*len*q6*sin(q2); 0 0];
        G = [0; mp*g*len*sin(q2)]; Bbar = [1 0; 0 0; 0 0; 0 1/m3];
        Mbar = [M zeros(2); zeros(2) eye(2)]; Cbar = [C zeros(2); zeros(2) zeros(2)];
        Gbar = [G; zeros(2,1)]; J = [-1 0; 0 0; 1/m2 -1/m2; 0 1/m3];
        Minv = inv(Mbar);
        cont = x(1:8,i); cont(4,1) = cont(4,1) - pi;
        u = KK*cont + LL*lambda;
        
        %Iterate the dynamics
        q_val(5:8,i+1) = q_val(5:8,i) + (-Minv*Gbar - Minv*Cbar*[q5;q6;q7;q8] + Minv*Bbar*u + Minv*J*lambda)*dt;
        q_val(1,i+1) = q_val(1,i) + q_val(5,i)*dt;
        q_val(2,i+1) = q_val(2,i) + q_val(6,i)*dt;
        q_val(3,i+1) = q_val(3,i) + q_val(7,i)*dt;
        q_val(4,i+1) = q_val(4,i) + q_val(8,i)*dt;
     
     %Coordinate changes
     x(1,i+1) = q_val(1,i+1); x(2,i+1) = q_val(3,i+1); x(3,i+1) = q_val(4,i+1); x(4,i+1) = q_val(2,i+1);
     x(5,i+1) = q_val(5,i+1); x(6,i+1) = q_val(7,i+1); x(7,i+1) = q_val(8,i+1); x(8,i+1) = q_val(6,i+1);
        
    end
    
    %CHECK STABILITY
    %check if it is unstable
    if max( abs(  x(:,end) ) ) >= 0.01
        counter2 = counter2 + 1;
    end
    
    %Print the figure if the contact force didnt exceed 10^7
    if flag == 0
        %hold on
        plot([0:dt:t+dt],x(1:4,:),'LineWidth',1)
        legend('x1','x2','x3','x4')
    end
    x_save = [x_save x(:,end)];
end

xh = x(1:8,1:end-1);

save('data.mat','lamh','xh')