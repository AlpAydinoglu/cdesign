%Calculates the gain matrices K, L and the Lyapunov function parameters P,
%Q, R for the cartpole model
%Make sure that PenBMI and Yalmip are in the MATLAB path

clear
clc
close all

%In this example, we solve the BMI feasibility problem starting from random
%initial points
counter = 0;
while counter<= 100
    clear
    counter = 0;
    counter = counter + 1;

%optimization parameters
eps = 10^-3; %epsilon
num_iter = 500; %number of iterations for the optimization problem
%bounds
bound_1 = 10; %bound on EBK
bound_2 = 10; %bound on EBL
bound_k = 50; %bound on K
bound_l = 20; %bound on L
bound = 0.1; %bound on the decay rate

%parameters
g = 9.81;
mp = 0.1;
mc = 1;
l = 0.5;
d1 = 0.1;
d2 = -0.1;

%problem data
%state matrix, input matrix, contact matrix
A =  [zeros(2) eye(2); 0 g*mp/mc 0 0; 0 g*(mc+mp)/(l*mc) 0 0]; B = [0; 0; 1/mc; 1/(l*mc)]; 
D = [zeros(3,2); 1/(l*mp) -1/(l*mp)];
Ec = [-1 l 0 0; 1 -l 0 0]; Fc = 0.1*eye(2); w = [d1;-d2];
%extract dimension information
n = size(A,2); %dimension of state space
k = size(B,2); %dimension of input
m = size(D,2); %number of contacts

%s-procedure matrices
%Ex + F \lambda >= 0 and \lambda >= 0
Z = [Ec Fc w; zeros(m,n) eye(m) zeros(m,1); zeros(1,m+n) 1];
%extension of Z to (x,\lambda,\bar{\lambda},\eta, 1)
Zjjj = [Ec Fc; zeros(m,n) eye(m)];
Zjj = [Zjjj zeros(2*m,2*m) [w;zeros(m,1)]];
Za = [Zjj; zeros(1,n+3*m) 1];
%\lambda_i E_i^T x + F_i^T \lambda = 0
for i = 1:m
    Zc_nonsym = [zeros(n+i-1, n+m+1); Ec(i,:) Fc(i,:) w(i,:); zeros(m-i+1, n+m+1)]; %create nonsymmetric version
    Zc{i} = (Zc_nonsym + Zc_nonsym')/2; %symmetrize the matrix
end
%\lambda_i E_i^T x + F_i^T \lambda = 0 (extension to (x,\lambda,\bar{\lambda},\eta,1)
for i = 1:m
    Zc_nonsym = [zeros(n+i-1, n+3*m+1); Ec(i,:) Fc(i,:) zeros(1,2*m) w(i,:); zeros(m-i, n+3*m+1); zeros(2*m+1, n+3*m+1)]; %create nonsymmetric version
    Zc_e{i} = (Zc_nonsym + Zc_nonsym')/2; %symmetrize the matrix
end
%\lambda^T * \eta = 0
for i = 1:m
   hold = zeros(n+3*m+1);
   hold(n+i,n+2*m+i) = 1;
   hold = (hold + hold')/2;
   Zgg{i} = hold;
end
%\kappa x^T x - z^T H z \geq 0
H = [Ec*A Ec*D Fc eye(m) zeros(m,1)];
Hu = H' * H;
O = sdpvar(m,n+3*m+1);
    
%variables
P = sdpvar(n,n);  Q = sdpvar(n,m); R = sdpvar(m,m);
W = sdpvar(2*m+1,2*m+1); U = sdpvar(2*m+1,2*m+1); 
tau_b = sdpvar(1,1); tau_h = sdpvar(1,1);
for i = 1:2*m
   tau{i} = sdpvar(1,1); 
end
for i=1:m
    tau_Zgg{i} = sdpvar(1,1);
end
K = sdpvar(k,n); L = sdpvar(k,m);

%constraints
N11 = A' * P + P * A + K' * B' * P + P * B * K;
N12 = P * B * L + P * D + A' * Q + K' * B' * Q;
N13 = Q;
N21 = N12';
N22 = L' * B' * Q + D' * Q + Q' * D  + Q' * B * L ;
N23 = R;
N31 = Q';
N32 = R;
N33 = zeros(m,m);
Nf = [N11 N12 N13; N21 N22 N23; N31 N32 N33];
N = [Nf zeros(n+2*m,m+1); zeros(m+1, n+3*m+1)];
M = [P Q zeros(n,1); Q' R zeros(m,1); zeros(1, n+m+1)];
%construct the inequalities
ineq1 = M - Z' * W * Z;
ineq2 = N + Za' * U * Za + O' * H/2 + H' * O/2;
for i = 1:m
    ineq1 = ineq1 - tau{i}*Zc{i};
    ineq2 = ineq2 + tau{m+i}*Zc_e{i} + tau_Zgg{i}*Zgg{i};
end
matl = eps*eye(n+m+1); matl(n+m+1,n+m+1) = 0;
%matk = zeros(n+3*m+1); matk(1:n,1:n)=bound*eye(n);
matk = [bound*P zeros(n,3*m+1); zeros(3*m+1, n+3*m+1)];
F = [ineq1 >= matl, ineq2 <= -matk];
%F = [F, ineq2 >= -0.5*eye(n+3*m+1)];
%F = [F, P <= 100*eye(4)];
F = [F, U(:) >= 0, W(:) >= 0, tau_h >= 0];
F = [F, [bound_k*eye(k) K; (K)' bound_k*eye(n)] >= 0];
F = [F, [bound_l*eye(k) L; (L)' bound_l*eye(m)] >= 0];

KKK = 10*(rand(1,4)-0.5);
assign(K,KKK); 
LLL = 10*(rand(1,2)-0.5);
assign(L,LLL); 

options=sdpsettings('solver', 'penbmi', 'penbmi.PBM_MAX_ITER', num_iter,'usex0',1,'savesolveroutput',1);
out_solver = optimize(F,[],options);

%stop if the desired feasibility criteria is met
feas = max ( out_solver.solveroutput.feas(3), out_solver.solveroutput.feas(2) );
if feas <= 10^-6
   break 
end
end

% display the resulting matrices
LL=double(L);
KK=double(K);

%save necessary matrices
save('controller.mat','KK','LL','A','B','D','n','m','k','Ec','Fc','w')