function dy = sys_affine(t,y,A,B,D,K,L,m,Fc,Ec,c,kappa,H,k)
tau = y(2);
yu = y(1);
lambda = pathlcp(Fc,Ec*yu+c+H*tau);
u = K*yu + L*lambda;
dtau = kappa*eye(k)*(u-tau);
dyu = A*yu + D*lambda + B*tau;
dy = [dyu;dtau];

if lambda' * lambda >= 1000
    lambda' * lambda;
end