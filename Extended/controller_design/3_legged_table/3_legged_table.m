clear all
clc
close all

%addpath PENBMI, YALMIP, MOSEK

%optimization parameters
eps = 10^-2;
num_iter = 100;

%parameters
m=1;
g=9.81;

mu1=0.1;
mu2=0.5;
mu3=1;

%problem data
%state matrix, input matrix, contact matrix
A = [0];
B = [1]/4; %input matrix
D = [0 1 -1 0 0 0]/4; %contact matrix
%complementarity constraints
Ec = [0;0;0;0;0;0];
Fc = [0 -1 -1 mu1 mu2 mu3; 1 1 -1 0 0 0; 1 -1 1 0 0 0; 0 0 0 1 1 1; 0 0 0 1 1 1; 0 0 0 1 1 1];
H = [0;1;-1; 0; 0; 0]; 
c = [0 ; 0; 0; -m*g; -m*g; -m*g];
kappa = 100;

%W SOL(q,F) is a singleton
W = [0 0 0 1 1 1];
p = size(W, 1);

%extract dimension information
n = size(A,2); %dimension of state space
k = size(B,2); %dimension of input
m = size(D,2); %number of contacts

%define variables
%the state variables
for i = 1:n
   x{i} = sdpvar(1,1); 
end
%variables related to the contact force
for i = 1:m
    lam{i} = sdpvar(1,1); %lambda (contact force)
    rho{i} = sdpvar(1,1); %rho (slack variable)
    lamd{i} = sdpvar(1,1); %\bar{\lambda} (directional der. of contact for.)
    tau{i} = sdpvar(1,1); %the low-pass filter dynamics
    mu{i} = sdpvar(1,1);
end

%basis vectors
%the state vector
x_basis = [];
for i = 1:n
    x_basis = [x_basis; x{i}];
end
%the contact vector and related variables
lam_basis = [];
lamd_basis = [];
tau_basis = [];
rho_basis = [];
mu_basis = [];
for i = 1:m
    lam_basis = [lam_basis; lam{i}];
    lamd_basis = [lamd_basis; lamd{i}];
    rho_basis = [rho_basis; rho{i}];
    mu_basis = [mu_basis; mu{i}];
end
for i = 1:k
   tau_basis = [tau_basis; tau{i}]; 
end
%general basis vector
basis = monolist([x_basis' lam_basis' lamd_basis' rho_basis' tau_basis' mu_basis'],1);
%size of the basis vector
sbasis = length(basis);

%Controllers
K = sdpvar(k,n,'full'); 
L = sdpvar(k,m,'full'); 

%Lyap Function
%Define the variables
P1 = sdpvar(n,n); P2 = sdpvar(n,p); P3 = sdpvar(p,p);
P4 = sdpvar(n,k); P5 = sdpvar(p,k); P6 = sdpvar(k,k);
p1 = sdpvar(1,n); p2 = sdpvar(1,m); p3 = sdpvar(1,k); z = sdpvar(1,1);
%Construct the Lyapunov function
V = x_basis' * P1 * x_basis + 2 * x_basis' * P2 * W * lam_basis ...
    + lam_basis' * W' * P3 * W * lam_basis + 2 * x_basis' * P4 * tau_basis ...
    + 2 * lam_basis' * W' * P5 * tau_basis + tau_basis' * P6 * tau_basis + z;
    %+ p1 * x_basis + p2 * lam_basis + p3 * tau_basis;
%Define the derivative vectors
xdot = A*x_basis + B*tau_basis + D*lam_basis;
taudot = kappa*eye(k)*K*x_basis + kappa*eye(k)*L*lam_basis - kappa*eye(k)*tau_basis;
%Define the derivative of the Lyapunov function
Vdot = 2 * x_basis' * P1 * xdot + 2* x_basis' * P2 * W * lamd_basis ...
    + 2 * lam_basis' * W' * P2' * xdot + 2 * lam_basis' * W' * P3 * W * lamd_basis ...
    + 2 * x_basis' * P4 * taudot + 2 * tau_basis' * P4' * xdot ...
    + 2 * lam_basis' * W' * P5 * taudot + 2 * tau_basis' * P5' * W * lamd_basis ...
    + 2 * tau_basis' * P6 * taudot;
    %+ p1 * xdot + p2 * lamd_basis + p3 * taudot;

%initalize the constraint set
F = [];

%S-procedure terms
%(Ex + F \lam + c + H \tau) >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    INEQ1{i} = Ec(i,:)*x_basis + Fc(i,:)*lam_basis + c(i) + H(i,:)*tau_basis; 
    INEQ1_M{i} = sdpvar(1,1); F = [F, INEQ1_M{i} >= 0];
    INEQ1_MD{i} = sdpvar(1,1); F = [F, INEQ1_MD{i} >= 0];
    INEQ1_MDD{i} = sdpvar(1,1); F = [F, INEQ1_MDD{i} >= 0];
end
%(\lam) >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    INEQ2{i} = lam_basis(i); INEQ2_M{i} = sdpvar(1,1); F = [F, INEQ2_M{i} >= 0];
    INEQ2_MD{i} = sdpvar(1,1); F = [F, INEQ2_MD{i} >= 0]; 
    INEQ2_MDD{i} = sdpvar(1,1); F = [F, INEQ2_MDD{i} >= 0];  
end
%(\lam)^2 >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
enum = 1:m;
cmb = combnk(enum, 2);
for i = 1:length(cmb)
        INEQ4{i} = lam_basis( cmb(i,1) )*lam_basis( cmb(i,2)  ); INEQ4_M{i} = sdpvar(1,1); F = [F, INEQ4_M{i} >= 0];
        INEQ4_M2{i} = sdpvar(1,1); F = [F, INEQ4_M2{i} >= 0]; 
        INEQ4_MM2{i} = sdpvar(1,1); F = [F, INEQ4_MM2{i} >= 0];
end

%(Ex + F \lam + c + H \tau)^2 >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(cmb)
        INEQ6{i} = (Ec(cmb(i,1),:)*x_basis + Fc(cmb(i,1),:)*lam_basis + c(cmb(i,1)) + H(cmb(i,1),:)*tau_basis)*(Ec(cmb(i,2) ,:)*x_basis + Fc(cmb(i,2) ,:)*lam_basis + c(cmb(i,2) ) + H(cmb(i,2) ,:)*tau_basis); 
        INEQ6_M{i} = sdpvar(1,1); F = [F, INEQ6_M{i} >= 0]; 
        INEQ6_M2{i} = sdpvar(1,1); F = [F, INEQ6_M2{i} >= 0];
        INEQ6_MM2{i} = sdpvar(1,1); F = [F, INEQ6_MM2{i} >= 0];
end

%\lam^T (Ex + F\lam + c + H\tau) = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    EQ1{i} = lam_basis(i) * INEQ1{i}; EQ1_M{i} = sdpvar(1,1); EQ1_MD{i} = sdpvar(1,1); EQ1_MDD{i} = sdpvar(1,1);
end

%E xdot + F \lamdot + H \taudot + rho = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    EQ2{i} = Ec(i,:)*xdot + Fc(i,:)*lamd_basis + H(i,:)*taudot + rho{i};
    a1{i} = sdpvar(1,n); a2{i} = sdpvar(1,m); a3{i} = sdpvar(1,k); 
    a4{i} = sdpvar(1,m); a5{i} = sdpvar(1,m); a6{i} = sdpvar(1,m);
    EQ2_MD{i} = a1{i}*x_basis + a2{i}*lam_basis + a3{i}*tau_basis...
        + a4{i}*rho_basis + a5{i}*lamd_basis + a6{i}*mu_basis;
end
%\lam_i \rho_i = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    EQ3{i} = lam{i} * rho{i}; EQ3_MD{i} = sdpvar(1,1);
end
%\bar{\lambda}_i + \mu_i = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    EQ4{i} = lamd_basis(i) + mu_basis(i);
    b1{i} = sdpvar(1,n); b2{i} = sdpvar(1,m); b3{i} = sdpvar(1,k); 
    b4{i} = sdpvar(1,m); b5{i} = sdpvar(1,m); b6{i} = sdpvar(1,m);
    EQ4_MD{i} = b1{i}*x_basis + b2{i}*lam_basis + b3{i}*tau_basis...
        + b4{i}*rho_basis + b5{i}*lamd_basis + b6{i}*mu_basis;
end
%(Ex + F \lambda + c + H \tau)_i \mu_i = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    EQ5{i} = (Ec(i,:)*x_basis + Fc(i,:)*lam_basis + c(i) + H(i,:)*tau_basis)*mu_basis(i); EQ5_MD{i} = sdpvar(1,1);
end
%\rho_i \mu_i = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    EQ6{i} = rho_basis(i)*mu_basis(i); EQ6_MD{i} = sdpvar(1,1);
end

% x^2 >= val
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INEQ_S = (x_basis' * x_basis) - 0.01; INEQ_S1 = sdpvar(1,1); INEQ_S2 = sdpvar(1,1);
F = [F , INEQ_S1 >= 0, INEQ_S2 >= 0];

alpha = sdpvar(1,1);
%inequalities
ineq1 = V - eps * (x_basis' * x_basis) - eps * (tau_basis' * tau_basis);
ineq2 = - Vdot; 
ineq3 = -V + alpha * (x_basis' * x_basis) + alpha * (tau_basis' * tau_basis);
%add S-procedure terms
for i = 1:m
    ineq1 = ineq1 - INEQ1{i}*INEQ1_M{i} - INEQ2{i}*INEQ2_M{i}  - EQ1{i}*EQ1_M{i};
    ineq2 = ineq2 - INEQ1{i}*INEQ1_MD{i} - INEQ2{i}*INEQ2_MD{i} ... 
        - EQ1{i}*EQ1_MD{i} - EQ2{i}*EQ2_MD{i} - EQ3{i}*EQ3_MD{i} ...
        - EQ4{i}*EQ4_MD{i} - EQ5{i}*EQ5_MD{i} - EQ6{i}*EQ6_MD{i};
    ineq3 = ineq3 - INEQ1{i}*INEQ1_MDD{i} - INEQ2{i}*INEQ2_MDD{i}  - EQ1{i}*EQ1_MDD{i};
end

for i = 1:length(cmb)
    ineq1 = ineq1 - INEQ4{i}*INEQ4_M{i} - INEQ6{i}*INEQ6_M{i};
    ineq2 = ineq2 - INEQ4{i}*INEQ4_M2{i} - INEQ6{i}*INEQ6_M2{i};
    ineq3 = ineq3 - INEQ4{i}*INEQ4_MM2{i} - INEQ6{i}*INEQ6_MM2{i};
end

F = [F, alpha >= 0.1];

%Construct the sos program
v = monolist([x_basis' lam_basis' tau_basis'],1);
Ke = sdpvar(length(v));
p_sos = v'*Ke*v;
F = [F, coefficients(ineq1-p_sos,[x_basis' lam_basis' tau_basis']) == 0, Ke>=0];

h = monolist([x_basis' lam_basis' lamd_basis' rho_basis' tau_basis' mu_basis'],1);
He = sdpvar(length(h));
q_sos = h'*He*h;
F = [F, coefficients(ineq2-q_sos,[x_basis' lam_basis' lamd_basis' rho_basis' tau_basis' mu_basis']) == 0, He>=0];

tp = monolist([x_basis' lam_basis' tau_basis'],1);
TPe = sdpvar(length(tp));
tp_sos = tp'*TPe*tp;
F = [F, coefficients(ineq3-tp_sos,[x_basis' lam_basis' tau_basis']) == 0, TPe>=0];

bound_l = 1;
F = [F, [bound_l*eye(k) L; (L)' bound_l*eye(m)] >= 0];
KKK = -200; 
assign(K,KKK);

%solve BMI
options=sdpsettings('solver', 'penbmi', 'penbmi.PBM_MAX_ITER', num_iter,'usex0',1,'savesolveroutput',1);
out_solver = optimize(F,[],options);

%CHECK SDP
% options = sdpsettings('solver','mosek');
% optimize(F, [], options);

% display the resulting matrices
LL=double(L);
KK=double(K);

%save the matrices
save('controller.mat','KK','LL','A','B','D','n','m','k','Ec','Fc','c','kappa','H')