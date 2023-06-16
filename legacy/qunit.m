function qunit(uf,Xn,Yn)

n2=round(length(uf)/2);
mag=sqrt(uf(1:n2).^2+uf(1+n2:end).^2);
magi=1./mag;
u1=[uf(1:n2).*magi;uf(1+n2:end).*magi];
x=reshape(Xn,[],1);
y=reshape(Yn,[],1);
quiver(x,y,u1(1:n2),u1(n2+1:end));
