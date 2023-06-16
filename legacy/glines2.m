function [ux,xx,yx,uy,xy,yy]=glines(un,xm,ym,xn,yn,jx,jy)

ux=jx*un;
uy=un*jy';

[xx,yx]=ndgrid(xm,yn);
[yy,xy]=ndgrid(ym,xn);

