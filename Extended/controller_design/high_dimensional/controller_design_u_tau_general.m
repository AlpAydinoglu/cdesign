clear all
clc
close all

%addpath
addpath(genpath('C:\Users\alp1a\OneDrive\Masaüstü\research\YALMIP-master\YALMIP-master'))
addpath 'C:\Users\alp1a\Dropbox\research\PenBMI\PENBMI2.1-LIMITED\matlab'
addpath 'C:\Program Files\Mosek\9.2\toolbox\R2015a'

% run = 0;
% 
% while( run <= 10)
%  
% clear all
% clc
% run = 0;
% run = run+1;

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
%K = [-95.2953  -90.8673  -63.7916  -18.7213    2.1605 -104.1293  -93.4882  -32.0898; -20.0027  -57.4949 -111.6282  -79.8556  -30.9723  -28.6045  -71.2010  -78.8693; -59.7504  -67.1223  -79.6496  -23.4354 -103.7441  -76.0156  -67.9578  -27.0616; -42.6220  -41.4139  -54.6301  -27.1910  -62.9535  -73.2837  -78.0740  -83.4867; -5.9657  -34.2569  -76.4876  -47.8869  -23.5592   -2.3569  -97.2464  -34.9039; -88.3858   -8.6324  -70.0914  -59.9688  -19.3174  -32.7515  -32.0486  -91.0160];
%L = [ 0.0020   -0.0176   -0.0494    0.0168    0.4074    0.0183    0.0412    0.3036   -0.0426   -0.2532; 0.0355   -0.0043   -0.3855   -0.1168    0.2621    0.0340    0.0137    0.2439   -0.0516    0.0428; 0.0077    0.0062    0.2203    0.0440    0.0918    0.1359    0.1590   -0.1657   -0.1151   -0.2698; -0.0023    0.0190   -0.2033    0.0529   -0.3153   -0.2657   -0.4297    0.1207    0.1315    0.4440; -0.0108   -0.0113    0.1111   -0.0136   -0.2309    0.0831    0.0854   -0.2478    0.1603    0.2660; -0.0300    0.0030    0.3128   -0.0049   -0.2789    0.0985    0.1657   -0.3764   -0.0498   -0.0937];

%Lyap Function
%Define the variables
P1 = sdpvar(n,n); P2 = sdpvar(n,p); P3 = sdpvar(p,p);
p1 = sdpvar(1,n); p2 = sdpvar(1,m);
%Construct the Lyapunov function
V = x_basis' * P1 * x_basis + 2 * x_basis' * P2 * W * lam_basis ...
    + lam_basis' * W' * P3 * W * lam_basis + p1 * x_basis; %+ p2 * lam_basis;
%Define the derivative vectors
u = K*x_basis + L*lam_basis; 
xdot = A*x_basis + B*u + D*lam_basis;
%Define the derivative of the Lyapunov function
Vdot = 2 * x_basis' * P1 * xdot + 2* x_basis' * P2 * W * lamd_basis ...
    + 2 * lam_basis' * W' * P2' * xdot + 2 * lam_basis' * W' * P3 * W * lamd_basis ...
    + p1 * xdot; %+ p2 * lamd_basis;

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

rng shuffle
%bound_k = 50;
bound_l = 1;
%F = [F, [bound_k*eye(k) K; (K)' bound_k*eye(n)] >= 0];
F = [F, [bound_l*eye(k) L; (L)' bound_l*eye(m)] >= 0];
%KKK = 1*randn(k,n);  
% KKK = -10*[1 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 0 1 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 1] - 100*rand(k,n);
%KKK = -50*rand(k,n);
%KKK = -105*rand(k,n);
%worked
LLL = 0*randn(k,m);
KKK = [-95.7801  -95.1695  -55.7481  -18.2896   -1.1267  -96.8082  -95.3561  -38.8248;
  -27.8772  -65.2360  -95.4123  -86.5410  -40.0750  -15.4781  -74.0096  -83.7490;
  -47.5990  -67.2228  -81.9349  -24.7823 -101.2143  -87.0657  -67.0684  -25.7182;
  -52.5962  -36.4966  -62.0396  -21.3055  -64.2054  -68.0781  -76.9591  -82.5445;
  -10.1578  -31.0900  -85.0588  -46.4695  -19.4662   -4.3753  -95.9772  -28.4519;
  -80.5702   -4.1201  -78.2399  -60.4749  -12.6558  -45.3071  -30.5034  -85.3218];

assign(K,KKK); assign(L,LLL);

%solve BMI
options=sdpsettings('solver', 'penbmi', 'penbmi.PBM_MAX_ITER', num_iter,'usex0',1,'savesolveroutput',1);
out_solver = optimize(F,[],options);
% feas = max ( out_solver.solveroutput.feas(3), out_solver.solveroutput.feas(2) );
% if feas <= 10^-6
%     break 
% end
% end

%options = sdpsettings('solver','penbmi', 'penbmi.PBM_MAX_ITER', num_iter);
% options = sdpsettings('solver','mosek');
% optimize(F, [], options);

% display the resulting matrices
LL=double(L);
KK=double(K);

%save the matrices
save('controller.mat','KK','LL','A','B','D','n','m','k','Ec','Fc','c')