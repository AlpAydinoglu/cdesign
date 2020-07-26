function dy = sys_affine(t,y,A,B,D,K,L,m,Fc,Ec,c,k)
lambda = pathlcp(Fc,Ec*y+c);
u = K*y + L*lambda;
dy = A*y + D*lambda + B*u;