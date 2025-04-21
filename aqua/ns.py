import matplotlib.pyplot as plt
import torch as pt
import numpy as np
import h5py

from .vis import vis
from .op2d import Op
from .problem import problem
from .advance import BDFEXT, Solution

from concurrent.futures import ThreadPoolExecutor
executor = ThreadPoolExecutor(max_workers=1)

def ns(Nx,Ny,tf,dt,cname,nu,u0=None,pfx="",ifconv=True,iotime=0,iptime=1,dtype=pt.float64,device=pt.device('cpu')):
    # usage: ns(32,32,100,1e-3,"ldc",1/1000,[]);

    bc="dddd";

    (ufun,f,dom,ifexact)=problem(cname);
    pd=problem(cname);
    ufun = pd['ufun']
    f = pd['f']
    dom=pd['dom']
    ifexact=pd['exact']

#   dtype=pt.float32
#   device=pt.device('mps')
    
    op = Op(Nx,Ny,dom,bc)
    Xn, Yn = np.meshgrid(op.x[0],op.x[1])
    uf = ufun(Xn,Yn,0,nu)
    uf = pt.tensor(uf,device=device,dtype=dtype)
    
    op.to(dtype=dtype,device=device)
    if u0 is not None:
        u0 = u0.to(device=device,dtype=dtype)
#       u = u0.u.to(device=device,dtype=dtype)
        uf = op.mask(u0.u) + op.maskc(uf)
        p = u0.p.to(device=device,dtype=dtype)
        if u0.T is not None:
            T = u0.T.to(device=device,dtype=dtype)
        else:
            T = None
        ic = u0
    else:
        ic = Solution(uf)
    
    iostep = int(np.round(iotime/dt))
    ipstep = int(np.round(iptime/dt))
    
    def cb(istep,t,sol,op):
        future = executor.submit(io_tasks,istep,t,sol,op)
        return future
#       io_tasks(istep,t,sol,op)
        
    def io_tasks(istep,t,sol,op):
        folder=pfx
        if np.mod(istep,ipstep) == 0:
            print(f'\nStep {istep}, t={t:.3f}')
            i=int(np.round((istep)/ipstep))
            vis(sol,op,i=i,folder=folder)
            
        if np.mod(istep,iostep) == 0:
            uc = sol.u.cpu()
            pc = sol.p.cpu()
            tc = t.cpu()
                
            i=int(np.round((istep)/iostep))
            istr = f'{i}'.rjust(6,'0')
            fname=f'{folder}snap_{istr}.h5'
            with h5py.File(fname,'w') as h5f:
                h5f.create_dataset('u',data=uc)
                h5f.create_dataset('p',data=pc)
                h5f.create_dataset('nu',data=nu)
                h5f.create_dataset('t',data=tc)
                if sol.T is not None:
                    h5f.create_dataset('T',data=sol.T.cpu())

    bdfext = BDFEXT(ic,op,tf,dt,nu)
    bdfext.to(device=device,dtype=dtype)

    sol = bdfext.forward(k=1,cb=cb,iostep=iostep,ipstep=ipstep)

    print('Done!')
    return sol,op
