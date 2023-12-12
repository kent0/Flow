function ub=pod(us,Bb);

u0=mean(us,2);

us0=us-u0;

ndump=size(us0,2);

G=zeros(ndump,ndump);

for j=1:ndump
    disp(['gengram ',num2str(j)]);
    Buj=Bb(us0(:,j));
    for i=1:ndump
        G(i,j)=dot(us0(:,i),Buj);
    end
end

G=.5*(G+G');

[V,D]=eig(G);

D=flip(diag(D));

V=flip(V,2);

ub=[u0 us0*V(:,1:end-1)];

for i=2:ndump
    ub(:,i)=ub(:,i)./sqrt(dot(ub(:,i),Bb(ub(:,i))));
end
