function dy = sys_affine(t,y,A,B,D,K,L,m,Fc,Ec,c,kappa,H,k)
tau = y(end-k+1:end);
yu = y(1:end-k);
lambda = pathlcp(Fc,Ec*yu+c+H*tau);
u = K*yu + L*lambda;
dtau = kappa*eye(k)*(u-tau);
dyu = A*yu + D*lambda + B*tau;
dy = [dyu;dtau];

if lambda' * lambda >= 1000
    lambda' * lambda;
end

if sum( abs(tau) ) <= 10^-3
   %test 
end