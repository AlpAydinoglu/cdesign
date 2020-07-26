clear
clc
close all

%addpath to use PenBMI
addpath 'C:\Users\Alp\Desktop\PenBMI\PENBMI2.1-LIMITED\matlab'

%optimization parameters
eps = 10^-3;
num_iter = 1000;
bound_1 = 10; %bound on EBK
bound_2 = 10; %bound on EBL
bound_x = 50; %bound on ||x||_2
bound = eps;

%parameters
g = 9.81;
mp = 0.5;
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

%calculated using calculate_lipschitz_affine
c = 136.4812;
d = 12.2041;

%set gamma and kappa
kappa_sqrt = bound_1 * bound_x + bound_2 * sqrt(c*bound_x^2 + d);
kappa = kappa_sqrt^2;
kappa = 1;

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
Hcons = zeros(n+3*m+1); Hcons(end,end) = kappa;
Huu = Hcons - Hu;

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
ineq2 = N + Za' * U * Za + tau_h*Huu;
for i = 1:m
    ineq1 = ineq1 - tau{i}*Zc{i};
    ineq2 = ineq2 + tau{m+i}*Zc_e{i} + tau_Zgg{i}*Zgg{i};
end
matl = eps*eye(n+m+1); matl(n+m+1,n+m+1) = 0;
matk = zeros(n+3*m+1); matk(1:n,1:n)=bound*eye(n);
F = [ineq1 >= matl, ineq2 <= -matk];
F = [F, U(:) >= 0, W(:) >= 0, tau_h >= 0];

% solve BMI
options=sdpsettings('solver', 'penbmi', 'penbmi.PBM_MAX_ITER', num_iter);
optimize(F,[],options);

% display the resulting matrices
LL=double(L);
KK=double(K);

%save necessary matrices
save('controller.mat','KK','LL','A','B','D','n','m','k','Ec','Fc','w')