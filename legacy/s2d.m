su=s2d(u,Bbx,Bby,Dx,Dy,Rx,Ry,ifupd);
persistent Sx Sy Dinv

if (ifupd)
    Ax=Dx'*Bhx*Dx; Ax=.5*(Ax+Ax');
    Ay=Dy'*Bhy*Dy; Ay=.5*(Ay+Ay');
    [Sx,lamx]=eig(Ax,Rx*Bbx*Rx');
    [Sy,lamy]=eig(Ay,Ry*Bby*Ry');
    Dinv=inv(kron(Ry*Ry',sparse(lamx))+kron(sparse(lamy),Rx*Rx'));
end

t1=a2u(By,Bx,a2u(Jpy,Jpx,p));
t2=[a2u(Iy,Dx',t1);a2u(Dy',Ix,t1)];

t3=a2u(Sy,Sx,Dinv,*a2u(Sy',Sx'));

su=a2u(Jpy',Jpyx',a2u(By,Bx,a2u(Iy,Dx,t3)+a2u(Dy,Ix,t3)));
