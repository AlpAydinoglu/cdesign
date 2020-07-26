function dy = sys_affine(t,y,A,B,D,K,L,G,J,m,Fc,Ec,w)
lambda = pathlcp(Fc,Ec*y+w);
dy = (A+B*K)*y + (D+B*L)*lambda;
for i = 1:m
   dy = dy + (G{i} + B*J{i})*y*lambda(i);
end