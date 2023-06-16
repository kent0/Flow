function uv=start_uv(x,y,t,n)

ns=1:n;
lam=(ns-.5)*pi;
T=exp(-lam.^2*t);
a=(48*(-1).^ns)./(pi*(2*ns-1)).^3;
u=cos(reshape(y,[],1)*lam)*(T.*a)'+1.5*(1-reshape(y,[],1).^2);
uv=[u;u*0];
