from ast import Call

from torch._prims_common import prod
from op2d import Op
from setops import setops
from problem import problem
from advance import BDFEXT
import matplotlib.pyplot as plt
import torch as pt
from util import mag, Solution
from vis import vis
import numpy as np
import h5py
import scipy.sparse.linalg as spla
from advance import Solution

def stokes(Nx,Ny,cname,u0=[],pfx="",heat=False):
    # usage: ns(32,32,100,1e-3,"ldc",1/1000,[]);

    bc="dddd";

#   u0,_,_,_=stokes(Nx,Ny,cname,pfx);

    rtol = 1e-6

    (ufun,f,dom,ifexact)=problem(cname);
    pd=problem(cname);
    ufun = pd['ufun']
    f = pd['f']
    dom=pd['dom']
    ifexact=pd['exact']

    dtype=pt.float64
    device=pt.device('cpu') # mps doesn't support fp64
    
    op = Op(Nx,Ny,dom,bc)

    nu=1
#   Xn = op.Xn.detach().cpu().numpy()
#   Yn = op.Yn.detach().cpu().numpy()
    Xn, Yn = np.meshgrid(op.x[0],op.x[1])
    uf = ufun(Xn,Yn,0,nu)
    uf = pt.tensor(uf)
    ub = op.maskc(uf)
    ru = -op.Ainv(op.R(op.Ab(ub)))
    rp = -op.D(ub+op.RT(ru))
    rp = rp.numpy().flatten()
    
    op.to(dtype=dtype,device=device)
    
    S = spla.LinearOperator((len(rp),len(rp)),matvec=op.S)
    
    class Counter():
        def __init__(self):
            self.i = 0
        def __call__(self,res=None):
            self.i+=1
            if res and np.mod(self.i+1,100) == 0:
                self.res = res
                print(f'Uzawa GMRES {self.i+1}: {res:.3e}')
        
    
    cb = Counter()
    p,info = spla.gmres(S,rp,rtol=rtol,callback=cb,callback_type='pr_norm',maxiter=8192)
    
    print(f'\nUzawa GMRES converged in {cb.i} iterations')
    
    Npx = op.Jxpn.shape[1]
    Npy = op.Jypn.shape[1]
    
    p = pt.tensor(p.reshape((Npy,Npx)))
    
    u = op.RT(ru + op.Ainv(op.DT(p,restrict=True)))
    u = u + ub
    
    # check residual
    ru = op.R(op.Ab(u) - op.DT(p))
    print(f'Velocity residual {pt.amax(ru):.3e}')
    rp = op.D(u)
    print(f'Pressure residual {pt.amax(rp):.3e}')

    print('Done!\n')
    
    if heat:
        def fop(x: np.ndarray):
            shape = x.shape
            y = x * 1.0
            a0 = pt.tensor(y)
            a1 = a0.reshape((op.Ry.shape[0],op.Rx.shape[0]))
            b = op.RT(a1)
            c = 0.01*op.Ab(b) + op.Cb(u,b)
            return op.R(c).reshape(shape).numpy()
            
        Ac_shape = (op.Rx.shape[0]*op.Ry.shape[0],op.Rx.shape[0]*op.Ry.shape[0])
        
        AC = spla.LinearOperator(Ac_shape,matvec=fop)
        
        cb = Counter()
        T = op.maskc(mag(u))
        r = -op.R(op.Ab(T))
        T0 = op.RT(op.Ainv(r))
        T += T0
            
#       Tb = op.maskc(T)
            
#       r = -op.R(0.01*op.Ab(Tb) + op.Cb(u,Tb))
#       r = r.flatten().numpy()
        
#       T,info = spla.gmres(AC,r,rtol=rtol,callback=cb,callback_type='pr_norm',maxiter=8192)
#       T = pt.tensor(T)
#       T = T.reshape((op.Ry.shape[0],op.Rx.shape[0]))
#       T = Tb + op.RT(T)
    else:
        T=None
        
    sol = Solution(u,p,T)
    return sol,op