function uv=st2_uv(x,y,t,w,nu);

xp=reshape(x,[],1);
yp=reshape(y,[],1);

u0=exp(-sqrt(w*.5/nu)*yp).*cos(w*t-sqrt(w*.5/nu)*yp);
v0=xp*0;

uv=[u0;v0];
