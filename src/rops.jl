function rops(Abx,Aby,Bbx,Bby,Dbx,Dby,Rx,Ry)
    Ax=Rx*Abx*Rx'; Ax=.5*(Ax+Ax');
    Ay=Ry*Aby*Ry'; Ay=.5*(Ay+Ay');

    Bx=Rx*Bbx*Rx'; Bx=.5*(Bx+Bx');
    By=Ry*Bby*Ry'; By=.5*(By+By');

    Dx=Rx*Dbx*Rx';
    Dy=Ry*Dby*Ry';

    Ix=speye(size(Rx,1));
    Iy=speye(size(Ry,1));

    return Ax,Ay,Bx,By,Dx,Dy,Ix,Iy
end
