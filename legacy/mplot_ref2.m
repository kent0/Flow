function mplot_ref2(un,xn,yn);

%close all

nx=length(xn);
ny=length(yn);

dx=(xn(nx)-xn(1))*nx/(nx-1);
x0=xn(1);

vn=reshape(un,nx,ny);

mx=nx*10;
my=ny*10;

xm=(zwf(mx)+1)*dx*.5+x0;
ym=zwuni(my-1);

jx=interp_f_mat(mx,nx);
jy=interp_mat(ym,yn);

um=jx*vn*jy';

s=surf(xm,ym,um','edgecolor','none');
hold on;

[ux,xx,yx,uy,xy,yy]=glines2(vn,xm,ym,xn,yn,jx,jy);

plot3(xx,yx,ux,'k');
plot3(xy,yy,uy,'k');
