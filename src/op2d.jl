function a2u(ay,ax,u)

#   fast-tensor implementation of (ay (x) ax) u

    mx,nx=size(ax); my,ny=size(ay);
    n=nx*ny;
    m=mx*my;
    dim=Int(round(size(u,1)/n));
    col=size(u,2);
    
    if col > 1
        au=zeros(m,dim,col);
        v=reshape(u,n,dim,col);
        for ic=1:col
        for id=1:dim
#           au[(1+(id-1)*m):(id*m),ic]=reshape(ax*reshape(u[(1+(id-1)*n):(id*n),ic],nx,ny)*ay',m,1);
            au[:,id,ic]=reshape(ax*reshape(v[:,id,ic],nx,ny)*ay',:);
        end
        end
        au=reshape(au,m*dim,col);
    elseif dim > 1
        au=zeros(m,dim);
        v=reshape(u,n,dim);
        for id=1:dim
            au[:,id]=reshape(ax*reshape(v[:,id],nx,ny)*ay',:);
        end
        au=reshape(au,:);
    else
        au=reshape(ax*reshape(u,nx,ny)*ay',:);
    end

    if isdiag(ax) && isdiag(ay)
        axd=Vector(diag(ax));
        ayd=Vector(diag(ay))';

        if col > 1
            au=zeros(m,dim,col);
            v=reshape(u,n,dim,col);
            for ic=1:col
            for id=1:dim
    #           au[(1+(id-1)*m):(id*m),ic]=reshape(ax*reshape(u[(1+(id-1)*n):(id*n),ic],nx,ny)*ay',m,1);
#               au[:,id,ic]=reshape(ax*reshape(v[:,id,ic],nx,ny)*ay',:);
                au[:,id,ic]=reshape(axd.*reshape(v[:,id,ic],nx,ny).*ayd,:);
            end
            end
            au=reshape(au,m*dim,col);
        elseif dim > 1
            au=zeros(m,dim);
            v=reshape(u,n,dim);
            for id=1:dim
                au[:,id]=reshape(axd.*reshape(v[:,id],nx,ny).*ayd,:);
            end
            au=reshape(au,:);
        else
            au=reshape(axd.*reshape(u,nx,ny).*ayd,:);
        end
    end

    return au;
end

function c2d(cn,uvn,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jx,Jy,Fx,Fy,Rx,Ry)
    udim=2;
    uvm=a2u(Jy,Jx,uvn);
    cm=a2u(Jy,Jx,cn);
    n=Int(round(length(uvm)/udim));

    ux=a2u(Iby,Dbx,uvm);
    uy=a2u(Dby,Ibx,uvm);

    if udim == 2
        t=[cm[1:n].*ux[1:n]+cm[1+n:end].*uy[1:n];
           cm[1:n].*ux[1+n:end]+cm[1+n:end].*uy[1+n:end]];
    else
        t=[cm[1:n].*ux[1:n]+cm[1+n:end].*uy[1:n]];
    end

    return a2u(Jy',Jx',a2u(Bby,Bbx,t));
end

function h2d(ub,Abx,Aby,Bbx,Bby,nu,bdti)
    Hx=nu*Abx+(bdti*.5)*Bbx
    Hy=nu*Aby+(bdti*.5)*Bby
    return a2u(Bby,Hx,ub)+a2u(Hy,Bbx,ub)
end

function dt2d(p,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,ifr)
    t=a2u(Bby,Bbx,a2u(Jpy,Jpx,p));
#   t=a2u(Jpy,Jpx,p);
    dtp=[a2u(Iby,Dbx',t);a2u(Dby',Ibx,t)];
    if ifr; dtp=a2u(Ry,Rx,dtp); end

    return dtp
end

let Sx,Sy,Dinv
    global function h2dinv(u,Ax,Ay,Bx,By,Ix,Iy,nu,bdti,ifupd)
        if ifupd
            Hx=nu*Ax+bdti*Bx*.5; Hx=.5*(Hx+Hx');
            Hy=nu*Ay+bdti*By*.5; Hy=.5*(Hy+Hy');
    #       Sx,lamx=eigen(Hx,full(Bx));
    #       Sy,lamy=eigen(Hy,full(By));
            lamx,Sx=eigen(Hx,Matrix(Bx));
            lamy,Sy=eigen(Hy,Matrix(By));
    #       Dinv=diag(inv(kron(Iy,sparse(lamx))+kron(sparse(lamy),Ix)));
            Dinv=Vector(diag(kinv(kron(Iy,spdiagm(lamx))+kron(spdiagm(lamy),Ix))));
        end
        mx,nx=size(Ax); my,ny=size(Ay);
        dim=Int(round(length(u)/(nx*ny)));
        t1=a2u(Sy',Sx',u);

        if (dim==2)
            n2=Int(round(length(u)/2));
            t2=[Dinv.*t1[1:n2];Dinv.*t1[n2+1:end]];
        else
            t2=Dinv.*t1;
        end

        return a2u(Sy,Sx,t2)
    end
end


function d2d(uv,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,ifrt)
    if ifrt; rtuv=a2u(Ry',Rx',uv); else; rtuv=uv; end
    n2=Int(round(length(rtuv)/2));
    u=rtuv[1:n2]; v=rtuv[n2+1:end];
#   du=a2u(Jpy',Jpx',a2u(Bby,Bbx,(a2u(Iby,Dbx,u)+a2u(Dby,Ibx,v))));
    return a2u(Jpy'*Bby,Jpx'*Bbx*Dbx,u)+a2u(Jpy'*Bby*Dby,Jpx'*Bbx,v);

#function du=d2d(uv,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,ifrt)
#    if (ifrt) rtuv=a2u(Ry',Rx',uv); else rtuv=uv; end
#    n2=round(length(rtuv)/2);
#    u=rtuv(1:n2); v=rtuv(n2+1:end);
#    du=a2u(Jpy',Jpx',a2u(Bby,Bbx,(a2u(Iby,Dbx,u)+a2u(Dby,Ibx,v))));
#end

end

let Sx=[],Sy=[],Dinv=[]
    global function sb2dinv(u,Bbx,Bby,Dbx,Dby,Jpx,Jpy,Rx,Ry);
    #persistent Dinv Sx Sy

        if size(Sx,1) != size(Rx,1)
            Bxi=kinv(Rx*Bbx*Rx');
            Byi=kinv(Ry*Bby*Ry');

            Bx=Jpx'*Bbx*Rx'*Bxi*Rx*Bbx*Jpx; Bx=.5*(Bx+Bx');
            By=Jpy'*Bby*Ry'*Byi*Ry*Bby*Jpy; By=.5*(By+By');

            Ax=Jpx'*Bbx*Dbx*Rx'*Bxi*Rx*Dbx'*Bbx*Jpx; Ax=.5*(Ax+Ax');
            Ay=Jpy'*Bby*Dby*Ry'*Byi*Ry*Dby'*Bby*Jpy; Ay=.5*(Ay+Ay');

            Ix=speye(size(Jpx,2));
            Iy=speye(size(Jpy,2));

            lamx,Sx=eigen(Ax,Bx);
            lamy,Sy=eigen(Ay,By);

            Dinv=1.0./diag(kron(spdiagm(lamy),Ix)+kron(Iy,spdiagm(lamx)));
        end

    #   tmp1=a2u(Sy',Sx',u);
    #   tmp2=Dinv.*tmp1;

        return a2u(Sy,Sx,Dinv.*a2u(Sy',Sx',u));
    end
end
