function siu=sb2dinv(u,Bbx,Bby,Dbx,Dby,Jpx,Jpy,Rx,Ry);
persistent Dinv Sx Sy

if size(Sx,1) ~= size(Rx,1)
    Bxi=inv(Rx*Bbx*Rx');
    Byi=inv(Ry*Bby*Ry');

    Bx=Jpx'*Bbx*Rx'*Bxi*Rx*Bbx*Jpx; Bx=.5*(Bx+Bx');
    By=Jpy'*Bby*Ry'*Byi*Ry*Bby*Jpy; By=.5*(By+By');

    Ax=Jpx'*Bbx*Dbx*Rx'*Bxi*Rx*Dbx'*Bbx*Jpx; Ax=.5*(Ax+Ax');
    Ay=Jpy'*Bby*Dby*Ry'*Byi*Ry*Dby'*Bby*Jpy; Ay=.5*(Ay+Ay');

    Ix=speye(size(Jpx,2));
    Iy=speye(size(Jpy,2));

    [Sx,lamx]=eig(Ax,Bx);
    [Sy,lamy]=eig(Ay,By);

    Dinv=1./diag(kron(lamy,Ix)+kron(Iy,lamx));
end

siu=a2u(Sy,Sx,Dinv.*a2u(Sy',Sx',u));
