%Calculates the gain matrices K, L and the Lyapunov function parameters P,
%Q, R for the 3box model
%Make sure that PenBMI and Yalmip are in the MATLAB path

clear
clc
close all

%optimization parameters
eps = 10^-3;
num_iter = 300;

%parameters
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

%s-procedure matrices
%Ex + F \lambda >= 0 and \lambda >= 0
Z = [Ec Fc; zeros(m,n) eye(m)];
%extension of Z to (x,\lambda,\bar{\lambda},\eta)
Za = [Z zeros(2*m,2*m)];
%\lambda_i E_i^T x + F_i^T \lambda = 0
for i = 1:m
    Zc_nonsym = [zeros(n+i-1, n+m); Ec(i,:) Fc(i,:); zeros(m-i, n+m)]; %create nonsymmetric version
    Zc{i} = (Zc_nonsym + Zc_nonsym')/2; %symmetrize the matrix
end
%\lambda_i E_i^T x + F_i^T \lambda = 0 (extension to (x,\lambda,\bar{\lambda},\eta)
for i = 1:m
    Zc_e{i} = [Zc{i} zeros(n+m,2*m); zeros(2*m,n+3*m)];
end
%\lambda^T * \eta = 0
for i = 1:m
   hold = zeros(n+3*m);
   hold(n+i,n+2*m+i) = 1;
   hold = (hold + hold')/2;
   Zgg{i} = hold;
end
%EA x + ED \lam + Fc \bar{\lam} + \eta = -EBK x - EBL \lam = 0
H = [Ec*A Ec*D Fc eye(m)];

%variables
P = sdpvar(n,n);  Q = sdpvar(n,m); R = sdpvar(m,m);
W = sdpvar(2*m,2*m); U = sdpvar(2*m,2*m); 
tau_b = sdpvar(1,1); O = sdpvar(m,n+3*m);
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
N = [Nf zeros(n+2*m,m); zeros(m, n+3*m)];
M = [P Q; Q' R];
%construct the inequalities
ineq1 = M - Z' * W * Z;
ineq2 = N + Za' * U * Za + O' * H/2 + H' * O / 2;
for i = 1:m
    ineq1 = ineq1 - tau{i}*Zc{i};
    ineq2 = ineq2 + tau{m+i}*Zc_e{i} + tau_Zgg{i}*Zgg{i};
end
matl = eye(n+m); matl(1:n,1:n) = eye(n); matk = zeros(n+3*m); matk(1:n,1:n)= eye(n);
F = [ineq1 >= eps*matl, ineq2 <= -0.01*matk];
F = [F, U(:) >= 0, W(:) >= 0];
F = [F, K(:,2) == 0, K(:,6) == 0]; %sparsity, no feedback on \ddot{x_2} and \dot{x_2}

% solve BMI
options=sdpsettings('solver', 'penbmi', 'penbmi.PBM_MAX_ITER', num_iter);
optimize(F,[],options);

% display the resulting matrices
LL=double(L);
KK=double(K);

%save necessary matrices
save('controller.mat','KK','LL','A','B','D','n','m','k','Ec','Fc')