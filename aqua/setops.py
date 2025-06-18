import torch as pt
from torch import Tensor
import math
from torch.nn import Parameter, Buffer

from .util import cf64

args = {
    'device': pt.device('cpu'), 
    'dtype': pt.float64,
    'requires_grad': False
}

def genops(N,dom,fullmass=False):
    ax=dom[0]
    bx=dom[1]
    ay=dom[2]
    by=dom[3]
    
    lx=bx-ax
    ly=by-ay
    
    Nx,Ny=N
    
    zxn,wxn=zwgll(Nx)
    zyn,wyn=zwgll(Ny)
    
    wxn=lx*wxn*.5; wyn=ly*wyn*.5;
    xn=ax+lx*(zxn+1)*.5; yn=ay+ly*(zyn+1)*.5;
    X=pt.stack((xn,yn), dim=0)
    
    if fullmass:
        zxn1,wxn1=zwgll(Nx+1)
        zyn1,wyn1=zwgll(Ny+1)
        Jxn1=J(zxn1,zxn)
        Jyn1=J(zyn1,zyn)
        
        Mbx1=pt.diag(wxn1)
        Mby1=pt.diag(wyn1)
            
        Mbx = Jxn1.T @ Mbx1 @ Jxn1
        Mby = Jyn1.T @ Mby1 @ Jyn1
    else:
        Mbx=wxn
        Mby=wyn
    
    Dbx=D(xn)
    Dby=D(yn)
    
    return Mbx,Mby,Dbx,Dby,X

def setops(Nu,dom,bc:str,fullmass=False):
    Np=(Nu[0]-2,Nu[1]-2)
    Nd=(int(math.ceil(1.5*Nu[0])),int(math.ceil(1.5*Nu[1])))
    
    Mbxn,Mbyn,Dbxn,Dbyn,Xn = genops(Nu,dom,fullmass=fullmass)
    Mbxm,Mbym,Dbxm,Dbym,Xm = genops(Nd,dom,fullmass=fullmass)
    Mbxp,Mbyp,Dbxp,Dbyp,Xp = genops(Np,dom,fullmass=fullmass)
    
    i1,i2,i3,i4 = 0,Nu[0],0,Nu[1]
    
    if bc[0] == 'd':
        i1+=1
        
    if bc[1] == 'd':
        i2-=1 
        
    if bc[2] == 'd':
        i3+=1
        
    if bc[3] == 'd':
        i4-=1
    
    Rx=pt.eye(Nu[0], **args)
    Rx=Rx[i1:i2,:]
    
    Ry=pt.eye(Nu[1], **args)
    Ry=Ry[i3:i4,:]
    
    if bc[1] == 'p' and bc[2] == 'p': Rx[0,-1]=1.
    if bc[3] == 'p' and bc[4] == 'p': Ry[0,-1]=1.
    
    Jxnm = J(Xm[0],Xn[0])
    Jynm = J(Xm[1],Xn[1])
    Jxpn = J(Xn[1],Xp[0])
    Jypn = J(Xn[1],Xp[1])
    
    return [
            Buffer(x) for x in [
            Mbxn,Mbyn,Dbxn,Dbyn,Xn,
            Mbxm,Mbym,Dbxm,Dbym,Xm,
            Mbxp,Mbyp,Dbxp,Dbyp,Xp,
            Jxnm,Jynm,Jxpn,Jypn,
            Rx,Ry,
        ]
    ]

def zwgc(n):
    """
    Computes the p+1 Gauss-Chebyshev nodes z on [-1,1]
    and the p+1 weights w on the cpu with float64's, gradient-free.
    """

    device = pt.device('cpu')
    dtype  = pt.float64

    k = pt.linspace(0, n-1, n, device=device, dtype=dtype)

    return -pt.cos((k+0.5) * (pt.pi / n))

def zwgll(n):
    """
    Computes the p+1 Gauss-Lobatto-Legendre nodes z on [-1,1]
    and the p+1 weights w on the cpu with float64's, gradient-free.
    """
    
    p = n - 1
    z = pt.zeros(n, **args)
    w = pt.zeros(n, **args)

    z[0] = -1
    z[-1] = 1

    if p > 1:
        if p == 2:
            z[1] = 0
        else:
            M = pt.zeros((p-1, p-1), **args)
            for i in range(p-2):
                M[i+1, i] = (i+1) * (i+3) / ((i+1.5) * (i+2.5))

            M = .5*pt.sqrt(M)
            
            D, V = pt.linalg.eigh(M)
            z[1:p],_ = D.sort()

    # Compute the weights w
    w[0] = 2 / (p * n)
    w[-1] = w[0]

    for i in range(1, p):
        x = z[i]
        z0 = 1
        z1 = x
        z2 = 0 # not necessary, but proper
        for j in range(1, p):
            z2 = x * z1 * (2 * j + 1) / (j + 1) - z0 * j / (j + 1)
            z0 = z1
            z1 = z2
        w[i] = 2 / (p * n * z2**2)

    return z, w

def J(xo: Tensor, xi: Tensor):
    xo = cf64(xo)
    xi = cf64(xi)
    
    no = len(xo)
    ni = len(xi)

    a = pt.ones(ni,**args)
    
    for i in range(ni):
        for j in range(i):
            a[i] *= (xi[i] - xi[j])
        for j in range(i + 1, ni):
            a[i] *= (xi[i] - xi[j])
    a = 1.0 / a

    J = pt.zeros((no, ni), **args)
    s = pt.ones(ni, **args)
    t = pt.ones(ni, **args)

    for i in range(no):
        x = xo[i]

        ind = pt.where(xi == x)[0]
        if len(ind) > 0:
            J[i, ind] = 1
            continue

        for j in range(1, ni):
            s[j] = s[j - 1] * (x - xi[j - 1])
            t[ni - j - 1] = t[ni - j] * (x - xi[ni - j])
        J[i, :] = a * s * t

    return J

def D(x):
    x = cf64(x)
    n=len(x)
    a=pt.ones(n, **args)
    d=pt.zeros(n,n, **args)
    
    for j in range(n):
        for i in range(n):
            d[i,j]=x[i]-x[j]
        d[j,j]=1
        
    
#   a[i] = pt.prod(d[i,:])
#   for i in range(n):
#       for j in range(i): 
#           a[i]*=x[i]-x[j]
#       for j in range((i+1),n):
#           a[i]*=x[i]-x[j]
            
    a = pt.prod(d,1)
    a=1/a; # These are the alpha_i's

    d=1/d
    
    for i in range(n):
        d[i,i]=0
        d[i,i]=pt.sum(d[i,:])

    for j in range(n):
        for i in range(n):
            if i!=j:
                d[i,j] = a[j]/(a[i]*(x[i]-x[j]))

    return d

#def legendre(x,N)
##   Compute the Legendre polynomials up to degree N evaluated at points x
#
#    m=len(x);
#    l=pt.ones(m,N+1) # zeroth order
#    l[:,1]=x         # first order
#
#    for n in range(1,N+1):
#        n = 0
#        i = n + 1
#        l[:,i] = ( (2*n+1)*x*l[:,i-1]-n*l[:,i-2] )/(n+1);
#
#    return l
#end
#
#
#function filt(z,wz,nf,df)
#    N=length(z)-1;
#
#    Lz=legendre(z,N);
#    Bh=spdiagm(wz);
#
#    for i=1:(N+1); Lz[:,i]=Lz[:,i]*sqrt(1/(Lz[:,i]'*Bh*Lz[:,i])); end
#    Fx=Lz*spdiagm([ones(N+1-nf);1.0.-df./(nf:-1:1)])*Lz'*Bh;
#
#    return Fx
#end
#
#
#function speye(N)
#    return sparse(I,N,N)
#end
#
#
#function ndgrid(xx,yy)
#    ox=ones(length(xx))
#    oy=ones(length(yy))
#
#    return kron(xx,oy'),kron(ox,yy')
#end
#
#
#function kinv(a)
#    if isdiag(a)
#        return spdiagm(1.0./Vector(diag(a)))
#    else
#        return inv(a)
#    end
#end
#
#function setplots#(Xn,Yn,Xn2,Yn2,Xm,Ym,Xpn,Ypn,Jxnl,Jynl,Jxn2l,Jyn2l,Jxml,Jyml,Jpxl,Jpyl,Xl,Yl,pdim,vort#,sf,Db)
#    n2=Int(round(length(Xn)));
#
#    plot_m1=(u,fname="figure")->mplot_ref(u,Xn,Yn,Jxnl,Jynl,Xl,Yl,pdim,fname);
#    plot_mm=(u,fname="figure")->mplot_ref(u,Xm,Ym,Jxml,Jyml,Xl,Yl,pdim,fname);
#    plot_m2=(u,fname="figure")->mplot_ref(u,Xpn,Ypn,Jpxl,Jpyl,Xl,Yl,pdim,fname);
#    plot_m3=(u,fname="figure")->mplot_ref(u,Xn2,Yn2,Jxn2l,Jyn2l,Xl,Yl,pdim,fname);
#
#    p1=(uv,fname="figure")->plot_m1(uv[1:n2],fname);
#    p2=(uv,fname="figure")->plot_m1(uv[n2+1:end],fname);
#    pm=(uv,fname="figure")->plot_m1(sqrt.(uv[1:n2].^2+uv[n2+1:end].^2),fname);
#    pv=(uv,fname="figure")->plot_m1(vort(uv),fname);
#    psf=(uv,ubsf,fname="figure")->plot_m1(sf(uv,ubsf),fname);
#
#    pp=(p,fname="figure")   ->plot_m2(p,fname);
#    pdiv=(uv,fname="figure")->plot_m2(Db(uv),fname);
#
#    return p1,p2,pm,pv,psf,pp,pdiv,plot_m1,plot_m2,plot_mm,plot_m3
#end

#function heaviside(x)
#    return 0.5*(sign.(x).+1);
#end
#
#
#function lag(alag,a)
#    res=alag
#    for i=0:size(res,2)-2
#        res[:,end-i]=res[:,end-i-1]
#    end
#    res[:,1]=a[:]
#
#    return res
#end
#
#
#function lag!(res,a);
#    for i=0:size(res,2)-2
#        res[:,end-i]=res[:,end-i-1]
#    end
#    res[:,1]=a
#end

#function evalf(uf,C,Hb,Binv,Einv,Db,DT,R,RT)
#    ub=uf-RT(R(uf));
#
#    f1=R(incomp2(RT(Binv(R(Hb(uf,0)+C(uf,uf)))),Binv,Einv,Db,DT,RT));
#    f2=Binv(DT(Einv(Db(ub))));
#
#    return f1,f2
#end
#
#function gcfl(uf,dt,Xn,Yn);
#    cflx=0.;
#    cfly=0.;
#    cflxy=0.;
#    n2=Int(round(length(uf)/2));
#    Nx,Ny=size(Xn); Nx=Nx-1; Ny=Ny-1;
#    Uf=reshape(uf[1:n2],Nx+1,Ny+1);
#    Vf=reshape(uf[n2+1:end],Nx+1,Ny+1);
#    for j=2:Ny
#    for i=2:Nx
#        minx=min(Xn[i+1,j]-Xn[i,j],Xn[i,j]-Xn[i-1,j]);
#        miny=min(Yn[i,j+1]-Yn[i,j],Yn[i,j]-Yn[i,j-1]);
#
#        cflxi=abs(Uf[i,j])*dt/minx;
#        cflyi=abs(Vf[i,j])*dt/miny;
#        cflxy=max(cflxy,cflxi+cflyi);
#    end
#    end
#    return cflxy
#end
#
#
#let z = 0
#    global function incrementZ()
#        z += 1
#        return z
#    end
#end
#
#let z = []
#    global function augmentZ()
#        append!(z,length(z))
#        return z
#    end
#end
#
#
#function interp_uv(uv,Nx,Ny,Mx,My)
#    zxn,_=zwgll(Nx);
#    zyn,_=zwgll(Ny);
#    zxm,_=zwgll(Mx);
#    zym,_=zwgll(My);
#    Jx=interp_mat(zxm,zxn);
#    Jy=interp_mat(zym,zyn);
#    Juv=a2u(Jy,Jx,uv);
#
#    return Juv
#end
