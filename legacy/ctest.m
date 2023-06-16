m=20;

%mop=@(s) R(Hb(RT(s),0)+C(uf,RT(s)));
mmat=amat+c2mat;

%A=kron(Ay,Bx)+kron(By,Ax);
xe=rand(size(mmat,1),1);
b=mmat*xe;

%mmat=opmat(aop,length(b))

[Q,H]=arnoldi(@(x)mmat*x,b,m);

es=eig(H(1:m,1:m));
scatter(real(es),imag(es),'x')

hold on;

[a,bb,c,x0]=ellipse(es);

d=x0;
c=c*1.;

niter=100;
[x1,rs1]=cheby(@(s)mmat*s,b,c,d,niter);
%[x2,rs2] = cheby2(H, b, b*0, niter, d*2, 0);
