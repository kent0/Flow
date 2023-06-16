function e=piter(A,n);

b=rand(size(A));

for i=1:n-1
    b=A*b;
end

e=norm(A*b)/norm(b);
