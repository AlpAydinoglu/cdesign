function dy = sys_general(t,y,A,B,D,K,L,Fc,Ec)
lambda = pathlcp(Fc,Ec*y);
dy = (A+B*K)*y + (D+B*L)*lambda;