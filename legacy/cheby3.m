function x=cheby(H,b,c,d,niter)

r=b;
rs(1)=norm(r);
del=(1/d)*r;
x=0*b;
x=x+del;

alpha=2*d/(2*d^2-c^2);
beta=d*alpha-1;

for i=1:niter
    r=b-H(x);
    del=alpha*r+beta*del;
    x=x+del;

    alpha=(d-(c*.5)^2*alpha)^(-1);
    beta=d*alpha-1;
end
