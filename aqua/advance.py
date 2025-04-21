from torch import Tensor
from torch import nn 
from torch.nn import Buffer
from setops import J, D
from typing import Optional
import torch as pt
import numpy as np
from util import Solution, mag
from torch import _VF

#def TimeStepper():
class BDFEXT(nn.Module):
    def __init__(self, ic, op, tf, dt, nu=1):
        super().__init__()
        self.ic = ic
        self.op = op
        self.tf = tf
        self.dt = dt
        self.nu = nu
        
        u = self.ic.u
        self.ulag = Buffer(pt.stack([u,u,u]))
        if hasattr(self.ic,'p') and len(self.ic.p) > 0:
            p = self.ic.p
        else:
            p = self.op.zero('p')
        self.plag = Buffer(pt.stack([p,p,p]))
        
        c = self.op.Cb(u,u)
        self.clag = Buffer(pt.stack([c,c,c]))
        
        if hasattr(self.ic,'T') and len(self.ic.T) > 0:
            T = self.ic.T
        else:
            T = self.op.zero('T')
            
        self.Tlag = Buffer(pt.stack([T,T,T]))
    
    def forward(self, k=3, cb=None,iostep=1e9,ipstep=1e9):
        
        args = {
            'dtype': self.ic.u.dtype,
            'device': self.ic.u.device,
            'requires_grad': False
        }
        
        nsteps = round(self.tf/self.dt)
        self.dt = self.tf / nsteps
        u = self.ic.u
        p = self.ic.p
        T = self.ic.T
        
        ts = pt.linspace(0,self.tf,nsteps+1,**args)
        td1 = lambda x,y: pt.tensordot(x,y, dims=1)
        ub = self.op.maskc(u)
        
        Tb = None
        if T is not None:
            Tb = mag(ub)
        
        bdti = []
        bdti1 = []
        alphas = []
        alphas2 = []
        Hbub=[]
        HbTb=[]
        
        if cb is not None:
            uc = u.clone()
            pc = p.clone()
            Tc = T.clone() if T is not None else None
            sol = Solution(uc,pc,Tc)
            cb(0,ts[0].clone(),sol,self.op)
            
        for i, t in enumerate(ts[1:]):
            o = min(k,i+1) # temporal order
            if i < 4:
                t3 = ts[max(0,i-o-1):i+2].flip(dims=[0])
                bdti,alphas,alphas2=bdfext_var(t3,3,**args)
                bdti1=bdti[0]
                self.op.set_Hb([bdti1,self.nu])
                self.op.Hinv.to(device=self.ic.u.device,dtype=self.ic.u.dtype)
                Hbub=self.op.Hb(ub)
                HbTb=self.op.Hb(Tb)
    
            # set up rhs
            
#           pt.tensordot(bdti[1:],self.ulag, dims=1)
            
            b=-self.op.Mb(td1(bdti[1:],self.ulag))
            b=b-td1(alphas,self.clag)
            b=b-Hbub
            
            ps=td1(alphas2,self.plag)
            b=b+self.op.DT(ps)
#           b=self.op.R(b+self.op.DT(ps))
#           us=self.op.RT(self.op.Hinv(b))+ub
            us=self.op.RT_Hinv_R(b)+ub
            u,p=self.op.incomp(us,ps,bdti1)
            
            if T is None:
                c = self.op.Cb(u,u)
            else:
                Tstar = td1(alphas,self.Tlag).reshape(-1,T.shape[-2],T.shape[-1])
                cmb = pt.cat([u,Tstar],dim=0)
                tmp = self.op.Cb(u,cmb)
                c = tmp[:-1]
                b=-self.op.Mb(td1(bdti[1:],self.Tlag))
                b=b-tmp[-1]
                b=b-HbTb
#               T=self.op.RT(self.op.Hinv(self.op.R(b)))+Tb
                T=self.op.RT_Hinv_R(b)+Tb
            
            self.lag(u,p,T,c)
            
            if cb is not None and (np.mod(i+1,ipstep) == 0 or np.mod(i+1,iostep) == 0):
                uc = u.clone()
                pc = p.clone()
                Tc = T.clone() if T is not None else None
                sol = Solution(uc,pc,Tc)
                cb(i+1,t.clone(),sol,self.op)
            
        sol = Solution(u,p,T)
        return sol
            
    def lag(self,u,p,T,c):
        for i in range(2):
            self.plag[2-i] = self.plag[1-i]
            self.ulag[2-i] = self.ulag[1-i]
            self.clag[2-i] = self.clag[1-i]
        self.plag[0] = p
        self.ulag[0] = u
        self.clag[0] = c
        
        if T is not None:
#           pt.roll(self.Tlag,shifts=(1,),dims=(0,))
            self.Tlag[1:] = self.Tlag[0:-1].clone()
            self.Tlag[0] = T
            
def bdfext_var(ts,k,dtype=None,device=None,requires_grad=False):
    alphas=(J(ts[0:1],ts[1:]).T).reshape(-1)
    alphas2=(J(ts[0:1],ts[1:-1]).T).reshape(-1)
    betas=(D(ts))[0,:]
    
    betas = pt.cat((betas,pt.zeros(k+1-len(betas))))
    alphas = pt.cat((alphas,pt.zeros(k-len(alphas))))
    alphas2 = pt.cat((alphas2,pt.zeros(k-len(alphas2))))
    
    betas = betas.to(dtype=dtype,device=device)
    alphas = alphas.to(dtype=dtype,device=device)
    alphas2 = alphas2.to(dtype=dtype,device=device)
    
    return betas,alphas,alphas2