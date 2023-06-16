function hiu=h2dinv(u,Ax,Ay,Bx,By,Ix,Iy,nu,bdti,ifupd)
    persistent Sx Sy Dinv
    if (ifupd)
        Hx=nu*Ax+bdti*Bx*.5; Hx=.5*(Hx+Hx');
        Hy=nu*Ay+bdti*By*.5; Hy=.5*(Hy+Hy');
        [Sx,lamx]=eig(Hx,full(Bx));
        [Sy,lamy]=eig(Hy,full(By));
        Dinv=diag(inv(kron(Iy,sparse(lamx))+kron(sparse(lamy),Ix)));
    end
    [mx,nx]=size(Ax); [my,ny]=size(Ay);
    dim=round(length(u)/(nx*ny));
    t1=a2u(Sy',Sx',u);

    if (dim==2)
        n2=round(length(u)/2);
        t2=[Dinv.*t1(1:n2);Dinv.*t1(n2+1:end)];
    else
        t2=Dinv.*t1;
    end

    hiu=a2u(Sy,Sx,t2);
end
