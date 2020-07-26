function dy = sys_affine(t,y,A,B,D,K,L,Fc,Ec,w)
lambda = pathlcp(Fc,Ec*y+w);
dy = (A+B*K)*y + (D+B*L)*lambda;