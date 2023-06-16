function [x,rs]=cheby(H,b,c,d,niter)

rs=zeros(niter+1,1);
r=b;
rs(1)=norm(r);
del=(1/d)*r;
x=0*b;
x=x+del;

alpha=2*d/(2*d^2-c^2);
beta=d*alpha-1;

for i=1:niter
    r=b-H(x);
    rs(i+1)=norm(r);
    del=alpha*r+beta*del;
    x=x+del;

    alpha=(d-(c*.5)^2*alpha)^(-1);
    beta=d*alpha-1;
end
