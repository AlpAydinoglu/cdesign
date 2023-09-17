clear all
clc
close all

%addpath YALMIP, PENBMI, MOSEK

%optimization parameters
eps = 10^-3;
num_iter = 50;

%problem data
%state matrix, input matrix, contact matrix
A = zeros(8);
B = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 1 0 0 0 0; 0 0 0 0 0 0; zeros(4,2) eye(4) ];
D = [0 0 0 0 -1 0 0 0 0 -1; -1 0 0 0 0 0 1 0 0 0; 0 0 1 0 0 0 0 0 0 1; 0 -1 0 0 0 0 0 1 0 0; 0 0 0 0 0 -1 0 0 -1 0; 0 0 0 0 0 0 -1 0 0 0; 0 0 0 1 0 0 0 0 1 0; 0 0 0 0 0 0 0 -1 0 0];
%complementarity constraints
Ec = [0 -1 0 0 0 0 0 0; 0 0 0 -1 0 0 0 0; 0 0 -1 0 0 0 0 0; 0 0 0 0 0 0 -1 0; 1 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0; 0 1 0 0 0 -1 0 0; 0 0 0 1 0 0 0 -1; 0 0 0 0 -1 0 1 0; -1 0 1 0 0 0 0 0 ];
Fc = eye(10);
c = zeros(10,1);

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
p1 = sdpvar(1,n); p2 = sdpvar(1,m);
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
for i = 1:m
    for j = 1:m
        ind = j+(i-1)*m;
        INEQ4{ind} = lam_basis(i)*lam_basis(j); INEQ4_M{ind} = sdpvar(1,1); F = [F, INEQ4_M{ind} >= 0];
        INEQ4_M2{ind} = sdpvar(1,1); F = [F, INEQ4_M2{ind} >= 0];
    end
end 

%(Ex + F \lam + c)^2 >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    for j = 1:m
        ind = j+(i-1)*m;
        INEQ6{ind} = (Ec(i,:)*x_basis + Fc(i,:)*lam_basis + c(i))*(Ec(j,:)*x_basis + Fc(j,:)*lam_basis + c(j)); 
        INEQ6_M{ind} = sdpvar(1,1); F = [F, INEQ6_M{ind} >= 0]; 
        INEQ6_M2{ind} = sdpvar(1,1); F = [F, INEQ6_M2{ind} >= 0];
    end
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
    %- eps * (lam_basis' * (W' * W) * lam_basis); 
ineq2 = - Vdot - (eps) * (x_basis' * x_basis);
%add S-procedure terms
for i = 1:m
    ineq1 = ineq1 - INEQ1{i}*INEQ1_M{i} - INEQ2{i}*INEQ2_M{i}  - EQ1{i}*EQ1_M{i};
    ineq2 = ineq2 - INEQ1{i}*INEQ1_MD{i} - INEQ2{i}*INEQ2_MD{i} ... 
        - EQ1{i}*EQ1_MD{i} - EQ2{i}*EQ2_MD{i} - EQ3{i}*EQ3_MD{i} ...
        - EQ4{i}*EQ4_MD{i} - EQ5{i}*EQ5_MD{i} - EQ6{i}*EQ6_MD{i};
end

for i = 1:(m*m)
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

%BOUNDS ON K and L
%bound_k = 50;
bound_l = 1;
%F = [F, [bound_k*eye(k) K; (K)' bound_k*eye(n)] >= 0];
F = [F, [bound_l*eye(k) L; (L)' bound_l*eye(m)] >= 0];
%An initialization that worked
LLL = zeros(k,m);
KKK = [-95.7801  -95.1695  -55.7481  -18.2896   -1.1267  -96.8082  -95.3561  -38.8248;
  -27.8772  -65.2360  -95.4123  -86.5410  -40.0750  -15.4781  -74.0096  -83.7490;
  -47.5990  -67.2228  -81.9349  -24.7823  -99.2143  -87.0657  -67.0684  -25.7182;
  -52.5962  -36.4966  -62.0396  -21.3055  -64.2054  -68.0781  -76.9591  -82.5445;
  -10.1578  -31.0900  -85.0588  -46.4695  -19.4662   -4.3753  -95.9772  -28.4519;
  -80.5702   -4.1201  -78.2399  -60.4749  -12.6558  -45.3071  -30.5034  -85.3218];

assign(K,KKK); assign(L,LLL);

%solve BMI
options=sdpsettings('solver', 'penbmi', 'penbmi.PBM_MAX_ITER', num_iter,'usex0',1,'savesolveroutput',1);
out_solver = optimize(F,[],options);

% display the resulting matrices
LL=double(L);
KK=double(K);

%save the matrices
save('controller.mat','KK','LL','A','B','D','n','m','k','Ec','Fc','c')
