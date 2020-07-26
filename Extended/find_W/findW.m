clear all
clc
close all

%addpath to use PenBMI
addpath(genpath('C:\Users\alp1a\OneDrive\Masaüstü\research\YALMIP-master\YALMIP-master'))
addpath 'C:\Program Files\Mosek\9.2\toolbox\R2015a'

%optimization parameters
l = 1; %dimension of rectangular matrix l x m

%complementarity constraints
Fc = [0 -1 -1; 1 1 -1; 1 -1 1]; 

%extract dimension information
m = size(Fc,2); %number of contacts

%define variables
%variables related to the contact force
for i = 1:m
    lam1{i} = sdpvar(1,1); %lambda1 (contact force)
    lam2{i} = sdpvar(1,1); %lambda2 (contact force)
    q{i} = sdpvar(1,1);
    x{i} = sdpvar(1,1);
end

%basis vectors
%the contact vector and related variables
lam_basis1 = [];
lam_basis2 = [];
q_basis = [];
for i = 1:m
    lam_basis1 = [lam_basis1; lam1{i}];
    lam_basis2 = [lam_basis2; lam2{i}];
    q_basis = [q_basis; q{i}];
end

%basis
basis = monolist([lam_basis1' lam_basis2' q_basis'],1);
sbasis = length(basis);

%Subspace parameters
A = sdpvar(l,m); 

%initalize the constraint set
F = [];

%S-procedure terms
%(q + F \lambda) >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    INEQ1{i} = q_basis(i) + Fc(i,:)*lam_basis1; 
    INEQ1_Ms{i} = sdpvar(sbasis,sbasis);
    INEQ1_M{i} = basis' * INEQ1_Ms{i} * basis; F = [F, INEQ1_Ms{i} >= 0]; 
    INEQ1_M2s{i} = sdpvar(sbasis,sbasis); 
    INEQ1_M2{i} = basis' * INEQ1_M2s{i} * basis; F = [F, INEQ1_M2s{i} >= 0];  
    INEQ11{i} = q_basis(i) + Fc(i,:)*lam_basis2; 
    INEQ11_Ms{i} = sdpvar(sbasis,sbasis);
    INEQ11_M{i} = basis' * INEQ11_Ms{i} * basis; F = [F, INEQ11_Ms{i} >= 0]; 
    INEQ11_M2s{i} = sdpvar(sbasis,sbasis);
    INEQ11_M2{i} = basis' * INEQ11_M2s{i} * basis; F = [F, INEQ11_M2s{i} >= 0];
end
%(\lam) >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    INEQ2{i} = lam_basis1(i);
    INEQ2_Ms{i} = sdpvar(sbasis,sbasis);
    INEQ2_M{i} = basis' * INEQ2_Ms{i} * basis; F = [F, INEQ2_Ms{i} >= 0];
    INEQ2_M2s{i} = sdpvar(sbasis,sbasis);
    INEQ2_M2{i} =  basis' * INEQ2_M2s{i} * basis; F = [F, INEQ2_M2s{i} >= 0];
    INEQ22{i} = lam_basis2(i); 
    INEQ22_Ms{i} = sdpvar(sbasis,sbasis);
    INEQ22_M{i} = basis' * INEQ22_Ms{i} * basis; F = [F, INEQ22_Ms{i} >= 0];
    INEQ22_M2s{i} = sdpvar(sbasis,sbasis);
    INEQ22_M2{i} = basis' * INEQ22_M2s{i} * basis; F = [F, INEQ22_M2s{i} >= 0];
end
%(\lam)^2 >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    for j = 1:m
        ind = j+(i-1)*m;
        INEQ4{ind} = lam_basis1(i)*lam_basis1(j);
        INEQ4_Ms{ind} = sdpvar(sbasis,sbasis);
        INEQ4_M{ind} = basis' * INEQ4_Ms{ind} * basis; F = [F, INEQ4_Ms{ind} >= 0];
        INEQ4_M2s{ind} = sdpvar(sbasis,sbasis);
        INEQ4_M2{ind} = basis' * INEQ4_M2s{ind} * basis; F = [F, INEQ4_M2s{ind} >= 0];
        INEQ44{ind} = lam_basis2(i)*lam_basis2(j); 
        INEQ44_Ms{ind} = sdpvar(sbasis,sbasis);
        INEQ44_M{ind} = basis' * INEQ44_Ms{ind} * basis; F = [F, INEQ44_Ms{ind} >= 0];
        INEQ44_M2s{ind} = sdpvar(sbasis,sbasis);
        INEQ44_M2{ind} = basis' * INEQ44_M2s{ind} * basis; F = [F, INEQ44_M2s{ind} >= 0];
        INEQ444{ind} = lam_basis1(i)*lam_basis2(j); 
        INEQ444_Ms{ind} = sdpvar(sbasis,sbasis);
        INEQ444_M{ind} = basis' * INEQ444_Ms{ind} * basis; F = [F, INEQ444_Ms{ind} >= 0];
        INEQ444_M2s{ind} = sdpvar(sbasis,sbasis);
        INEQ444_M2{ind} = basis' * INEQ444_M2s{ind} * basis; F = [F, INEQ444_M2s{ind} >= 0];
    end
end

%(q + F \lambda)^2 >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    for j = 1:m
        ind = j+(i-1)*m;
        INEQ6{ind} = (q_basis(i) + Fc(i,:)*lam_basis1)*(q_basis(j) + Fc(j,:)*lam_basis1);
        INEQ6_Ms{ind} = sdpvar(sbasis,sbasis);
        INEQ6_M{ind} = basis' * INEQ6_Ms{ind} * basis; F = [F, INEQ6_Ms{ind} >= 0];
        INEQ6_M2s{ind} = sdpvar(sbasis,sbasis);
        INEQ6_M2{ind} = basis' * INEQ6_M2s{ind} * basis; F = [F, INEQ6_M2s{ind} >= 0];
        INEQ66{ind} = (q_basis(i) + Fc(i,:)*lam_basis2)*(q_basis(j) + Fc(j,:)*lam_basis2);
        INEQ66_Ms{ind} = sdpvar(sbasis,sbasis);
        INEQ66_M{ind} = basis' * INEQ66_Ms{ind} * basis; F = [F, INEQ66_Ms{ind} >= 0]; 
        INEQ66_M2s{ind} = sdpvar(sbasis,sbasis);
        INEQ66_M2{ind} = basis' * INEQ66_M2s{ind} * basis; F = [F, INEQ66_M2s{ind} >= 0];
        INEQ666{ind} = (q_basis(i) + Fc(i,:)*lam_basis1)*(q_basis(j) + Fc(j,:)*lam_basis2); 
        INEQ666_Ms{ind} = sdpvar(sbasis,sbasis);
        INEQ666_M{ind} = basis' * INEQ666_Ms{ind} * basis; F = [F, INEQ666_Ms{ind} >= 0]; 
        INEQ666_M2s{ind} = sdpvar(sbasis,sbasis);
        INEQ666_M2{ind} = basis' * INEQ666_M2s{ind} * basis; F = [F, INEQ666_M2s{ind} >= 0];
    end
end

%(q + F \lambda)_i lambda_j >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    for j = 1:m
        ind = j+(i-1)*m;
        INEQ9{ind} = (q_basis(i) + Fc(i,:)*lam_basis1)*lam_basis2(j); 
        INEQ9_Ms{ind} = sdpvar(sbasis,sbasis);
        INEQ9_M{ind} = basis' * INEQ9_Ms{ind} * basis; F = [F, INEQ9_Ms{ind} >= 0]; 
        INEQ9_M2s{ind} = sdpvar(sbasis,sbasis);
        INEQ9_M2{ind} = basis' * INEQ9_M2s{ind} * basis; F = [F, INEQ9_M2s{ind} >= 0];
        INEQ99{ind} = (q_basis(i) + Fc(i,:)*lam_basis2)*lam_basis1(j); 
        INEQ99_Ms{ind} = sdpvar(sbasis,sbasis);
        INEQ99_M{ind} = basis' * INEQ99_Ms{ind} * basis; F = [F, INEQ99_Ms{ind} >= 0]; 
        INEQ99_M2s{ind} = sdpvar(sbasis,sbasis);
        INEQ99_M2{ind} = basis' * INEQ9_M2s{ind} * basis; F = [F, INEQ99_M2s{ind} >= 0];
    end
end

%\lam^T (q + Fc \lam) = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    EQ1{i} = lam_basis1(i) * INEQ1{i};
    EQ1_Ms{i} = sdpvar(sbasis,sbasis); EQ1_M{i} = basis' * EQ1_Ms{i} * basis;
    EQ1_M2s{i} = sdpvar(sbasis,sbasis); EQ1_M2{i} = basis' * EQ1_M2s{i} * basis;
    EQ11{i} = lam_basis2(i) * INEQ11{i};
    EQ11_Ms{i} = sdpvar(sbasis,sbasis); EQ11_M{i} = basis' * EQ11_Ms{i} * basis;
    EQ11_M2s{i} = sdpvar(sbasis, sbasis); EQ11_M2{i} = basis' * EQ11_M2s{i} * basis;
end
    
%(q ) >= 0 (TRIAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    INEQ8{i} = q_basis(i)*q_basis(i) - 0.1; 
    INEQ8_Ms{i} = sdpvar(sbasis,sbasis);
    INEQ8_M{i} = basis' * INEQ8_Ms{i} * basis; F = [F, INEQ8_Ms{i} >= 0]; 
    INEQ8_M2s{i} = sdpvar(sbasis,sbasis); 
    INEQ8_M2{i} = basis' * INEQ8_M2s{i} * basis; F = [F, INEQ8_M2s{i} >= 0]; 
end
%INEQ8{1} = q_basis(1);
%INEQ8{2} = q_basis(2);
%INEQ8{3} = -q_basis(3);

pol_lam1 = lam_basis1' * lam_basis1 + lam_basis2' * lam_basis2;
pol_lam2 = lam_basis1' * lam_basis1 + lam_basis2' * lam_basis2;

epsilon = sdpvar(1,1);
%inequalities
ineq1 =  ( A*(lam_basis1 - lam_basis2 ) )*(pol_lam1) + epsilon * (pol_lam1);
ineq2 =  -( A*(lam_basis1 - lam_basis2 ) )*(pol_lam2) + epsilon * (pol_lam2);

%add S-procedure terms
for i = 1:m
    ineq1 = ineq1 - INEQ1{i}*INEQ1_M{i} - INEQ2{i}*INEQ2_M{i}  - EQ1{i}*EQ1_M{i}...
                  - INEQ11{i}*INEQ11_M{i} - INEQ22{i}*INEQ22_M{i} - EQ11{i}*EQ11_M{i};
    ineq2 = ineq2 - INEQ1{i}*INEQ1_M2{i} - INEQ2{i}*INEQ2_M2{i} - EQ1{i}*EQ1_M2{i}...
                  - INEQ11{i}*INEQ11_M2{i} - INEQ22{i}*INEQ22_M2{i} - EQ11{i}*EQ11_M2{i};
end

for i = 1:(m*m)
    ineq1 = ineq1 - INEQ4{i}*INEQ4_M{i} - INEQ6{i}*INEQ6_M{i}...
                  - INEQ44{i}*INEQ44_M{i} - INEQ66{i}*INEQ66_M{i}...
                  - INEQ444{i}*INEQ444_M{i} - INEQ666{i}*INEQ666_M{i}...
                  - INEQ9{i}*INEQ9_M{i} - INEQ99{i}*INEQ99_M{i};
    ineq2 = ineq2 - INEQ4{i}*INEQ4_M2{i} - INEQ6{i}*INEQ6_M2{i}...
                  - INEQ44{i}*INEQ44_M2{i} - INEQ66{i}*INEQ66_M2{i}...
                  - INEQ444{i}*INEQ444_M2{i} - INEQ666{i}*INEQ666_M2{i}...
                  - INEQ9{i}*INEQ9_M2{i} - INEQ99{i}*INEQ99_M2{i};
end

%Construct the sos program
v = monolist([lam_basis1' lam_basis2' q_basis'],2);
Ke = sdpvar(length(v));
p_sos = v'*Ke*v;
F = [F, coefficients(ineq1-p_sos,[lam_basis1' lam_basis2' q_basis']) == 0, Ke>=0];

h = monolist([lam_basis1' lam_basis2' q_basis'],2);
He = sdpvar(length(h));
q_sos = h'*He*h;
F = [F, coefficients(ineq2-q_sos,[lam_basis1' lam_basis2' q_basis']) == 0, He>=0];

F = [F, A(:) <= 1, A(:) >= -1, epsilon <= 0.001];

rng(1); %for consistency

%first step of the algorithm
r = rand(1,3);
obj = r*A';

%for the second step of algorithm
%r = rand(1,2);
%N = [-0.7071    0.7071; 0.5000    0.5000; 0.5000    0.5000];
%obj = r * N' * A' ;

%solve the optimization
%obj = [];
options = sdpsettings('solver', 'mosek');
mosek_opt = optimize(F,obj,options);

double(A)