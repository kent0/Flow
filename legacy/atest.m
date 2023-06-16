m=round(.05*size(jmat,1));
e1=eig(acmat);
r=rand(size(acmat,1),1);
[Q,H]=arnoldi(acmat,r,m); es=eig(H(1:m,1:m));
hold off; scatter(real(e1),imag(e1));
hold on; scatter(real(es),imag(es),'*');
