function aiu=a2dinv(u,Ax,Ay,Bx,By,Rx,Ry)
    persistent Sx Sy Dinv
    if length(Dinv) ~= length(Ax)
        [Sx,lamx]=eig(Ax,full(Bx));
        [Sy,lamy]=eig(Ay,full(By));
        Ix=Rx*Rx';
        Iy=Ry*Ry';
        Dinv=full(diag(inv(kron(Iy,sparse(lamx))+kron(sparse(lamy),Ix))));
    end
    [mx,nx]=size(Ax); [my,ny]=size(Ay);
    dim=round(length(u)/(nx*ny));
    t1=a2u(Sy',Sx',u);

    if (dim==2)
        t2=reshape(Dinv.*reshape(t1,nx*ny,2),[],1);
    else
        t2=Dinv.*t1;
    end

    aiu=a2u(Sy,Sx,t2);
end
