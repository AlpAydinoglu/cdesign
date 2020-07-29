function dy = sys_affine(t,y,A,B,D,K,L,m,Fc,Ec,c,kappa,H,k)
tau = y(2);
yu = y(1);
%lambda = pathlcp(Fc,Ec*yu+c+H*tau);

mu1=0.1;
mu2=0.5;
mu3=1;
if abs(yu) >= 0.7
    rng(1)
    lambda(4) = 9.81*rand(1);
    lambda(5) = (9.81-lambda(4))*rand(1);
    lambda(6) = 9.81-lambda(4)-lambda(5);
elseif abs(yu) < 0.7 && abs(yu) >= 0.5
    rng(2)
    lambda(4) = 9.81*rand(1);
    lambda(5) = (9.81-lambda(4))*rand(1);
    lambda(6) = 9.81-lambda(4)-lambda(5);
elseif abs(yu) < 0.5 && abs(yu) >= 0.3
    rng(3)
    lambda(4) = 9.81*rand(1);
    lambda(5) = (9.81-lambda(4))*rand(1);
    lambda(6) = 9.81-lambda(4)-lambda(5);
else
    rng(4)
    lambda(4) = 9.81*rand(1);
    lambda(5) = (9.81-lambda(4))*rand(1);
    lambda(6) = 9.81-lambda(4)-lambda(5);
end

lambda(1:3) = pathlcp([0 -1 -1; 1 1 -1; 1 -1 1], [mu1*lambda(4) + mu2*lambda(5) + mu3*lambda(6); tau; -tau]);
lambda = lambda';

u = K*yu + L*lambda;
dtau = kappa*eye(k)*(u-tau);
dyu = A*yu + D*lambda + B*tau;
dy = [dyu;dtau];