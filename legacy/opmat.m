function amat=opmat(aop,n)

z=zeros(n,1);
amat=zeros(n,n);

for i=1:n
    z(max(i-1,1))=0.;
    z(i)=1;
    amat(:,i)=aop(z);
end
