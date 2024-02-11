using SparseArrays

function zwgll(p);

#   computes the p+1 Gauss-Lobatto-Legendre nodes z on [-1,1]
#   i.e. the zeros of the first derivative of the Legendre polynomial
#   of degree p plus -1 and 1
#   and the p+1 weights w

    n=p+1;
    z=zeros(n);
    w=zeros(n);

    z[1]=-1;
    z[n]=+1;

    if p > 2
        M=zeros(p-1,p-1);
        for i=1:p-2
            M[i,i+1]=(1/2)*sqrt((i*(i+2))/((i+1/2)*(i+3/2)));
            M[i+1,i]=M[i,i+1];
        end
      
        D=eigvals(M);
        z[2:p]=sort(D);
    end

    # compute the weights w

    w[1]=2/(p*n);
    w[n]=w[1];

    for i=2:p
        x=z[i];

        z0=1;
        z1=x;
        z2=0;
        for j=1:p-1
            z2=x*z1*(2*j+1)/(j+1)-z0*j/(j+1);
            z0=z1;
            z1=z2;
        end;
        w[i]=2/(p*(n)*z2*z2);
    end

    return z,w

end


function interp_mat(xo,xi)

#   Compute the interpolation matrix from xi to xo

    no = length(xo);
    ni = length(xi);
    a  = ones(ni);
    for i=1:ni
        for j=1:(i-1);  a[i]=a[i]*(xi[i]-xi[j]); end
        for j=(i+1):ni; a[i]=a[i]*(xi[i]-xi[j]); end
    end
    a=1.0./a;

    J = zeros(no,ni);
    s = ones(ni);
    t = ones(ni);
    for i=1:no;
        x=xo[i];
        for j=2:ni
            s[j]      = s[j-1]*(x-xi[j-1]);
            t[ni+1-j] = t[ni+2-j]*(x-xi[ni+2-j])
        end
        J[i,:]=a.*s.*t;
    end

    return J
end


function deriv_mat(x)

#   Compute the interpolation matrix from x to x

    ni = length(x);
    a  = ones(ni);
    d  = zeros(ni,ni);
    for i=1:ni
        for j=1:(i-1);  a[i]=a[i]*(x[i]-x[j]); end
        for j=(i+1):ni; a[i]=a[i]*(x[i]-x[j]); end
    end
    a=1.0./a; # These are the alpha_i's

    for j=1:ni; for i=1:ni; d[i,j]=x[i]-x[j]; end; d[j,j]=1; end;
    d=1.0./d;
    for i=1:ni; d[i,i]=0; d[i,i]=sum(d[i,:]); end;

    for j=1:ni; for i=1:ni;
        if i!=j; d[i,j] = a[j]/( a[i]*(x[i]-x[j])); end;
    end;end;


    return d
end


function legendre(x,N)

#   Compute the Legendre polynomials up to degree N evaluated at points x

    m=length(x);
    l=ones(m,N+1);
    l[:,2]=x;

    for k=2:N
        i=k+1;
        l[:,i] = ( (2*k-1)*x.*l[:,i-1]-(k-1)*l[:,i-2] )/k;
    end

    return l
end


function filt(z,wz,nf,df)
    N=length(z)-1;

    Lz=legendre(z,N);
    Bh=spdiagm(wz);

    for i=1:(N+1); Lz[:,i]=Lz[:,i]*sqrt(1/(Lz[:,i]'*Bh*Lz[:,i])); end
    Fx=Lz*spdiagm([ones(N+1-nf);1.0.-df./(nf:-1:1)])*Lz'*Bh;

    return Fx
end


function speye(N)
    return sparse(I,N,N)
end


function ndgrid(xx,yy)
    ox=ones(length(xx))
    oy=ones(length(yy))

    return kron(xx,oy'),kron(ox,yy')
end


function kinv(a)
    if isdiag(a)
        return spdiagm(1.0./Vector(diag(a)))
    else
        return inv(a)
    end
end

function setplots(Xn,Yn,Xn2,Yn2,Xm,Ym,Xpn,Ypn,Jxnl,Jynl,Jxn2l,Jyn2l,Jxml,Jyml,Jpxl,Jpyl,Xl,Yl,pdim,vort,sf,Db)
    n2=Int(round(length(Xn)));

    plot_m1=(u)->mplot_ref(u,Xn,Yn,Jxnl,Jynl,Xl,Yl,pdim);
    plot_mm=(u)->mplot_ref(u,Xm,Ym,Jxml,Jyml,Xl,Yl,pdim);
    plot_m2=(u)->mplot_ref(u,Xpn,Ypn,Jpxl,Jpyl,Xl,Yl,pdim);
    plot_m3=(u)->mplot_ref(u,Xn2,Yn2,Jxn2l,Jyn2l,Xl,Yl,pdim);

    p1=(uv)->plot_m1(uv[1:n2]);
    p2=(uv)->plot_m1(uv[n2+1:end]);
    pm=(uv)->plot_m1(sqrt.(uv[1:n2].^2+uv[n2+1:end].^2));
    pv=(uv)->plot_m1(vort(uv));
    psf=(uv,ubsf)->plot_m1(sf(uv,ubsf));

    pp=(p)   ->plot_m2(p);
    pdiv=(uv)->plot_m2(Db(uv));

    return p1,p2,pm,pv,psf,pp,pdiv,plot_m1,plot_m2,plot_mm,plot_m3
end

function heaviside(x)
    return 0.5*(sign.(x).+1);
end


function lag(alag,a)
    res=alag
    for i=0:size(res,2)-2
        res[:,end-i]=res[:,end-i-1]
    end
    res[:,1]=a[:]

    return res
end


function lag!(res,a);
    for i=0:size(res,2)-2
        res[:,end-i]=res[:,end-i-1]
    end
    res[:,1]=a
end

function incomp(usf,ps,bdti1,Binv,Einv,Db,DT,RT)
    dp=-ortho(Einv(Db(usf))*bdti1);
    p=ps+dp;
    uf=usf+RT(Binv(DT(dp))*(1.0/bdti1));

    return uf,p
end

function incomp2(usf,Binv,Einv,Db,DT,RT)
    return usf+RT(Binv(DT(-ortho(Einv(Db(usf))))));
end

function ortho(p)
    return p.-mean(p)
end


function evalf(uf,C,Hb,Binv,Einv,Db,DT,R,RT)
    ub=uf-RT(R(uf));

    f1=R(incomp2(RT(Binv(R(Hb(uf,0)+C(uf,uf)))),Binv,Einv,Db,DT,RT));
    f2=Binv(DT(Einv(Db(ub))));

    return f1,f2
end

function gcfl(uf,dt,Xn,Yn);
    cflx=0.;
    cfly=0.;
    cflxy=0.;
    n2=Int(round(length(uf)/2));
    Nx,Ny=size(Xn); Nx=Nx-1; Ny=Ny-1;
    Uf=reshape(uf[1:n2],Nx+1,Ny+1);
    Vf=reshape(uf[n2+1:end],Nx+1,Ny+1);
    for j=2:Ny
    for i=2:Nx
        minx=min(Xn[i+1,j]-Xn[i,j],Xn[i,j]-Xn[i-1,j]);
        miny=min(Yn[i,j+1]-Yn[i,j],Yn[i,j]-Yn[i,j-1]);

        cflxi=abs(Uf[i,j])*dt/minx;
        cflyi=abs(Vf[i,j])*dt/miny;
        cflxy=max(cflxy,cflxi+cflyi);
    end
    end
    return cflxy
end


let z = 0
    global function incrementZ()
        z += 1
        return z
    end
end

let z = []
    global function augmentZ()
        append!(z,length(z))
        return z
    end
end


function interp_uv(uv,Nx,Ny,Mx,My)
    zxn,_=zwgll(Nx);
    zyn,_=zwgll(Ny);
    zxm,_=zwgll(Mx);
    zym,_=zwgll(My);
    Jx=interp_mat(zxm,zxn);
    Jy=interp_mat(zym,zyn);
    Juv=a2u(Jy,Jx,uv);

    return Juv
end

