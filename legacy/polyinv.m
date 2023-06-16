function b = polyinv(A,x,n);

c=x;
b=c*1;

for i=1:n
    c=c-A(c);
    b=b+c;
end
