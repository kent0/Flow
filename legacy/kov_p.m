function p=kov_uv(x,y,nu);

xp=reshape(x,[],1);
yp=reshape(y,[],1);

r=1./nu;
lam=.5*r-sqrt(r^2*.25+4*pi*pi);
u0=1-exp(lam*xp).*cos(2*pi*yp);
v0=(lam*.5/pi)*exp(lam*xp).*sin(2*pi*yp);

uv=[u0;v0];
p=-.5*exp(xp*(1-sqrt(1+16*pi*pi*nu*nu))/nu);
