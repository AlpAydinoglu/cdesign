clear all
clc
close all

%addpath to use PenBMI
addpath 'C:\Users\Alp\Desktop\PenBMI\PENBMI2.1-LIMITED\matlab'

%optimization parameters
eps = 10^-3;
num_iter = 100;

%problem data
%state matrix, input matrix, contact matrix
%A = [-0.1163   -0.6715; -0.3099   -0.4712];
A = rand(2);
B = [1 0;0 1]; %input matrix
%D = -eye(2); %contact matrix
D = eye(2);
%complementarity constraints
Ec = -[1 0; 0 1];
%Fc = [1 0; 0 2]; %Fc = eye(2);
%Fc = rand(2); 
Fc = eye(2);
%eig(Fc)
%Fc = [0.9595    0.0357; 0.6557    0.8491];
H = [1 0;0 0];
c = [1;1];
kappa = 100;
bound_k = 10;
bound_l = 10;

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
K = sdpvar(k,n,'full'); L = sdpvar(k,m,'full'); 

%Lyap Function
%Define the variables
P1 = sdpvar(n,n); P2 = sdpvar(n,m); P3 = sdpvar(m,m);
P4 = sdpvar(n,k); P5 = sdpvar(m,k); P6 = sdpvar(k,k);
p1 = sdpvar(1,n); p2 = sdpvar(1,m); p3 = sdpvar(1,k);
%Construct the Lyapunov function
V = x_basis' * P1 * x_basis + 2 * x_basis' * P2 * lam_basis ...
    + lam_basis' * P3 * lam_basis + 2 * x_basis' * P4 * tau_basis ...
    + 2 * lam_basis' * P5 * tau_basis + tau_basis' * P6 * tau_basis ...
    + p1 * x_basis + p2 * lam_basis + p3 * tau_basis;
%Define the derivative vectors
xdot = A*x_basis + B*tau_basis + D*lam_basis;
taudot = kappa*eye(k)*K*x_basis + kappa*eye(k)*L*lam_basis - kappa*eye(k)*tau_basis;
%Define the derivative of the Lyapunov function
Vdot = 2 * x_basis' * P1 * xdot + 2* x_basis' * P2 * lamd_basis ...
    + 2 * lam_basis' * P2' * xdot + 2 * lam_basis' * P3 * lamd_basis ...
    + 2 * x_basis' * P4 * taudot + 2 * tau_basis' * P4' * xdot ...
    + 2 * lam_basis' * P5 * taudot + 2 * tau_basis' * P5' * lamd_basis ...
    + 2 * tau_basis' * P6 * taudot + p1 * xdot + p2 * lamd_basis + p3 * taudot;

%initalize the constraint set
F = [];

%S-procedure terms
%(Ex + F \lam + c + H \tau) >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    INEQ1{i} = Ec(i,:)*x_basis + Fc(i,:)*lam_basis + c(i) + H(i,:)*tau_basis; 
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
    INEQ4{i} = lam_basis(i)^2; INEQ4_M{i} = sdpvar(1,1); F = [F, INEQ4_M{i} >= 0];
    INEQ4_MD{i} = sdpvar(1,1); F = [F, INEQ4_MD{i} >= 0];    
end
%%%%%WRITTEN IN A FAST MANNER FOR NOW
INEQ5 = lam_basis(1)*lam_basis(2); INEQ5_M = sdpvar(1,1); F = [F, INEQ5_M >= 0];
INEQ5_MD = sdpvar(1,1); F = [F, INEQ5_MD >= 0];  

%(Ex + F \lam + c + H \tau)^2 >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%in1 = H*K; in2 = H*L;
for i = 1:m
    INEQ6{i} = (Ec(i,:)*x_basis + Fc(i,:)*lam_basis + c(i) + H(i,:)*tau_basis)^2; 
    INEQ6_M{i} = sdpvar(1,1); F = [F, INEQ6_M{i} >= 0];
    INEQ6_MD{i} = sdpvar(1,1); F = [F, INEQ6_MD{i} >= 0];
end
%%%%%WRITTEN IN A FAST MANNER FOR NOW
INEQ7 = INEQ1{1}*INEQ1{2}; INEQ7_M = sdpvar(1,1); F = [F, INEQ7_M >= 0];
INEQ7_MD = sdpvar(1,1); F = [F, INEQ7_MD >= 0];

%\lam^T (Ex + F\lam + c + H\tau) = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    EQ1{i} = lam_basis(i) * INEQ1{i}; EQ1_M{i} = sdpvar(1,1); EQ1_MD{i} = sdpvar(1,1);
end
%E xdot + F \lamdot + H \taudot + rho = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    for i = 1:m
    EQ2{i} = Ec(i,:)*xdot + Fc(i,:)*lamd_basis + H(i,:)*taudot + rho{i};
    a1{i} = sdpvar(1,n); a2{i} = sdpvar(1,m); a3{i} = sdpvar(1,k); 
    a4{i} = sdpvar(1,m); a5{i} = sdpvar(1,m); a6{i} = sdpvar(1,m);
    EQ2_MD{i} = a1{i}*x_basis + a2{i}*lam_basis + a3{i}*tau_basis...
        + a4{i}*rho_basis + a5{i}*lamd_basis + a6{i}*mu_basis;
    end
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
    EQ5{i} = (Ec(i,:)*x_basis + Fc(i,:)*lam_basis + c(i) + H(i,:)*tau_basis)*mu_basis(i); 
    EQ5_MD{i} = sdpvar(1,1);
end
%\rho_i \mu_i = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    EQ6{i} = rho_basis(i)*mu_basis(i); EQ6_MD{i} = sdpvar(1,1);
end

%inequalities
ineq1 = V - eps * x_basis' * x_basis - eps * tau_basis' * tau_basis...
    - INEQ5*INEQ5_M - INEQ7*INEQ7_M - eps * lam_basis' * lam_basis;
ineq2 = - Vdot - INEQ5*INEQ5_MD - INEQ7*INEQ7_MD - eps * x_basis' * x_basis;

%add S-procedure terms
for i = 1:m
    ineq1 = ineq1 - INEQ1{i}*INEQ1_M{i} - INEQ2{i}*INEQ2_M{i}  - EQ1{i}*EQ1_M{i} ... 
                - INEQ4{i}*INEQ4_M{i} - INEQ6{i}*INEQ6_M{i};
    ineq2 = ineq2 - INEQ1{i}*INEQ1_MD{i} - INEQ2{i}*INEQ2_MD{i} ... 
        -INEQ4{i}*INEQ4_MD{i} - EQ1{i}*EQ1_MD{i} - EQ2{i}*EQ2_MD{i} - EQ3{i}*EQ3_MD{i} ...
        -INEQ6{i}*INEQ6_MD{i} - EQ4{i}*EQ4_MD{i} - EQ5{i}*EQ5_MD{i} - EQ6{i}*EQ6_MD{i};
end

%F = [F, EQ4_MD{2} >= 0.1];

%Construct the sos program
v = monolist([x_basis' lam_basis' tau_basis'],1);
Ke = sdpvar(length(v));
p_sos = v'*Ke*v;
F = [F, coefficients(ineq1-p_sos,[x_basis' lam_basis' tau_basis']) == 0, Ke>=0];

h = monolist([x_basis' lam_basis' lamd_basis' rho_basis' tau_basis' mu_basis'],1);
He = sdpvar(length(h));
q_sos = h'*He*h;
F = [F, coefficients(ineq2-q_sos,[x_basis' lam_basis' lamd_basis' rho_basis' tau_basis' mu_basis']) == 0, He>=0];

%solve BMI
options=sdpsettings('solver', 'penbmi', 'penbmi.PBM_MAX_ITER', num_iter,'usex0',1);
optimize(F,[],options);

% display the resulting matrices
LL=double(L);
KK=double(K);

%save the matrices
save('controller.mat','KK','LL','A','B','D','n','m','k','Ec','Fc','c','kappa','H')