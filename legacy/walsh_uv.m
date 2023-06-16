function uv=walsh_uv(x,y,t,nu);
% psi = (1/4) cos(3x) sin(4y) - (1/5) cos(5y) - (1/5) sin(5x)
% u = psi_y, v = -psi_x

% psi = cos(x) cos(y)

xp=reshape(x,[],1);
yp=reshape(y,[],1);

lam=-25;
u0=cos(3*xp).*cos(4*yp)+sin(5*yp);
v0=(3/4)*sin(3*xp).*sin(4*yp)+cos(5*xp);

%lam=-2;
%u0=-sin(xp).*cos(yp);
%v0=cos(xp).*sin(yp);

uv=[u0;v0]*exp(nu*lam*t);
