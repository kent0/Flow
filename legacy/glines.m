function [ux,xmn,ymn,uy,xnm,ynm]=glines(vn,xn,yn,Jx,Jy)

un=reshape(vn,[],1);

Ix=speye(size(Jx,2));
Iy=speye(size(Jy,2));

ux=a2u(Iy,Jx,un);
uy=a2u(Jy,Ix,un);

X=reshape(xn,[],1);
Y=reshape(yn,[],1);

xnm=a2u(Iy,Jx,X);
ynm=a2u(Iy,Jx,Y);

xmn=a2u(Jy,Ix,X);
ymn=a2u(Jy,Ix,Y);
