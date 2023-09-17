%This code is for systems with c >= 0, and a = 0 (we take z = 0, p = 0, r = 0). If
%one has a system with c < 0 or a \neq 0, add an upperbound on the Lyapunov
%function V similar to ineq1. Feel free to reach out to me if there are
%problems ( alpayd@seas.upenn.edu ).

clear all
clc
close all

%addpath
%here addpath of YALMIP, PENBMI and potentially MOSEK if you want to verify
%your controller using an SDP

%optimization parameters
eps = 10^-3;
num_iter = 50;

%problem data (fill in the matrices)
%state matrix, input matrix, contact matrix
A = [];
B = [];
D = [];
%complementarity constraints
Ec = [];
Fc = [];
c = [];

%extract dimension information
n = size(A,2); %dimension of state space
k = size(B,2); %dimension of input
m = size(D,2); %number of contacts

%W SOL(q,F) is a singleton
W = eye(m);
p = size(W, 1);

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
rho_basis = [];
mu_basis = [];
for i = 1:m
    lam_basis = [lam_basis; lam{i}];
    lamd_basis = [lamd_basis; lamd{i}];
    rho_basis = [rho_basis; rho{i}];
    mu_basis = [mu_basis; mu{i}];
end

%general basis vector
basis = monolist([x_basis' lam_basis' lamd_basis' rho_basis' mu_basis'],1);
%size of the basis vector
sbasis = length(basis);

%Controllers
K = sdpvar(k,n,'full'); 
L = sdpvar(k,m,'full');

%Lyap Function
%Define the variables
P1 = sdpvar(n,n); P2 = sdpvar(n,p); P3 = sdpvar(p,p);
%Construct the Lyapunov function
V = x_basis' * P1 * x_basis + 2 * x_basis' * P2 * W * lam_basis ...
    + lam_basis' * W' * P3 * W * lam_basis;
%Define the derivative vectors
u = K*x_basis + L*lam_basis; 
xdot = A*x_basis + B*u + D*lam_basis;
%Define the derivative of the Lyapunov function
Vdot = 2 * x_basis' * P1 * xdot + 2* x_basis' * P2 * W * lamd_basis ...
    + 2 * lam_basis' * W' * P2' * xdot + 2 * lam_basis' * W' * P3 * W * lamd_basis;

%initalize the constraint set
F = [];

%S-procedure terms
%(Ex + F \lam + c ) >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    INEQ1{i} = Ec(i,:)*x_basis + Fc(i,:)*lam_basis + c(i); 
    INEQ1_M{i} = sdpvar(1,1); F = [F, INEQ1_M{i} >= 0];
    INEQ1_MD{i} = sdpvar(1,1); F = [F, INEQ1_MD{i} >= 0];
end
%(\lam) >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    INEQ2{i} = lam_basis(i); INEQ2_M{i} = sdpvar(1,1); F = [F, INEQ2_M{i} >= 0];
    INEQ2_MD{i} = sdpvar(1,1); F = [F, INEQ2_MD{i} >= 0];    
end
%(\lam)^2 >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
enum = 1:m;
cmb = combnk(enum, 2);
for i = 1:length(cmb)
        INEQ4{i} = lam_basis( cmb(i,1) )*lam_basis( cmb(i,2)  ); INEQ4_M{i} = sdpvar(1,1); F = [F, INEQ4_M{i} >= 0];
        INEQ4_M2{i} = sdpvar(1,1); F = [F, INEQ4_M2{i} >= 0];
end

%(Ex + F \lam + c)^2 >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(cmb)
        INEQ6{i} = (Ec(cmb(i,1),:)*x_basis + Fc(cmb(i,1),:)*lam_basis + c(cmb(i,1)) )*(Ec(cmb(i,2) ,:)*x_basis + Fc(cmb(i,2) ,:)*lam_basis + c(cmb(i,2) ) ); 
        INEQ6_M{i} = sdpvar(1,1); F = [F, INEQ6_M{i} >= 0]; 
        INEQ6_M2{i} = sdpvar(1,1); F = [F, INEQ6_M2{i} >= 0];
end

%\lam^T (Ex + F\lam + c) = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    EQ1{i} = lam_basis(i) * INEQ1{i}; EQ1_M{i} = sdpvar(1,1); EQ1_MD{i} = sdpvar(1,1);
end

%E xdot + F \lamdot + rho = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    EQ2{i} = Ec(i,:)*xdot + Fc(i,:)*lamd_basis + rho{i};
    a1{i} = sdpvar(1,n); a2{i} = sdpvar(1,m); 
    a4{i} = sdpvar(1,m); a5{i} = sdpvar(1,m); a6{i} = sdpvar(1,m);
    EQ2_MD{i} = a1{i}*x_basis + a2{i}*lam_basis...
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
    b1{i} = sdpvar(1,n); b2{i} = sdpvar(1,m); 
    b4{i} = sdpvar(1,m); b5{i} = sdpvar(1,m); b6{i} = sdpvar(1,m);
    EQ4_MD{i} = b1{i}*x_basis + b2{i}*lam_basis...
        + b4{i}*rho_basis + b5{i}*lamd_basis + b6{i}*mu_basis;
end
%(Ex + F \lambda + c)_i \mu_i = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    EQ5{i} = (Ec(i,:)*x_basis + Fc(i,:)*lam_basis + c(i))*mu_basis(i); EQ5_MD{i} = sdpvar(1,1);
end
%\rho_i \mu_i = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    EQ6{i} = rho_basis(i)*mu_basis(i); EQ6_MD{i} = sdpvar(1,1);
end

%inequalities
ineq1 = V - eps * (x_basis' * x_basis);
ineq2 = - Vdot - (eps) * (x_basis' * x_basis); %comment out (eps) * (x_basis' * x_basis) for Lyapunov stability
%add S-procedure terms
for i = 1:m
    ineq1 = ineq1 - INEQ1{i}*INEQ1_M{i} - INEQ2{i}*INEQ2_M{i}  - EQ1{i}*EQ1_M{i};
    ineq2 = ineq2 - INEQ1{i}*INEQ1_MD{i} - INEQ2{i}*INEQ2_MD{i} ... 
        - EQ1{i}*EQ1_MD{i} - EQ2{i}*EQ2_MD{i} - EQ3{i}*EQ3_MD{i} ...
        - EQ4{i}*EQ4_MD{i} - EQ5{i}*EQ5_MD{i} - EQ6{i}*EQ6_MD{i};
end

for i = 1:length(cmb)
    ineq1 = ineq1 - INEQ4{i}*INEQ4_M{i} - INEQ6{i}*INEQ6_M{i};
    ineq2 = ineq2 - INEQ4{i}*INEQ4_M2{i} - INEQ6{i}*INEQ6_M2{i};
end

%Construct the sos program
v = monolist([x_basis' lam_basis'],1);
Ke = sdpvar(length(v));
p_sos = v'*Ke*v;
F = [F, coefficients(ineq1-p_sos,[x_basis' lam_basis']) == 0, Ke>=0];

h = monolist([x_basis' lam_basis' lamd_basis' rho_basis' mu_basis'],1);
He = sdpvar(length(h));
q_sos = h'*He*h;
F = [F, coefficients(ineq2-q_sos,[x_basis' lam_basis' lamd_basis' rho_basis' mu_basis']) == 0, He>=0];

%if bounds on ||K|| \leq a_1 and ||L|| \leq a_2 are preferred (details can
%be found on conference paper)
%bound_k = a_1;
%bound_l = a_2;
%F = [F, [bound_k*eye(k) K; (K)' bound_k*eye(n)] >= 0];
%F = [F, [bound_l*eye(k) L; (L)' bound_l*eye(m)] >= 0];

%initialization for matrices K and L (if commented, PenBMI chooses how to
%initialize them
% LLL = [];
% KKK = [];
% assign(K,KKK); assign(L,LLL);

%solve BMI
options=sdpsettings('solver', 'penbmi', 'penbmi.PBM_MAX_ITER', num_iter,'usex0',1,'savesolveroutput',1);
out_solver = optimize(F,[],options);

%If one wants to verify using an SDP solver after fixing the controllers K
%and L
% options = sdpsettings('solver','mosek');
% optimize(F, [], options);

% display the resulting matrices
LL=double(L);
KK=double(K);

%save the matrices
save('controller.mat','KK','LL','A','B','D','n','m','k','Ec','Fc','c')
