function S=orthoadd(S,s);

p=s;
for i=1:size(S,2);
    p=p-(p'*S(:,i))*S(:,i);
end
S=[S p/norm(p)];
