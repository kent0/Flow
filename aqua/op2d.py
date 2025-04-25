import torch as pt
from torch import Tensor, _VF
from torch.nn import Buffer
import scipy as sp
import numpy as np
import torch.backends.opt_einsum as opt_einsum
from typing import Any

from .setops import setops, J
from .util import cf64

def ortho(p):
    return p-pt.mean(p)
    
#opt_einsum.enabled = False
#pt.backends.
opt_einsum.strategy = 'optimal'

path_dict = {}

def einsum(*args: Any) -> Tensor:
    # taken from torch/functional.py
    import torch.backends.opt_einsum as opt_einsum

    # This wrapper exists to support variadic args.
    if len(args) < 2:
        raise ValueError(
            "einsum(): must specify the equation string and at least one operand, "
            "or at least one operand and its subscripts list"
        )

    equation = None
    operands = None

    if isinstance(args[0], Tensor):
        def parse_subscript(n: int) -> str:
            if n == Ellipsis:
                return "..."
            if n >= 0 and n < 26:
                return chr(ord("A") + n)
            if n >= 26 and n < 52:
                return chr(ord("a") + n - 26)
            raise ValueError(
                "einsum(): subscript in subscript list is not within the valid range [0, 52)"
            )

        equation = ",".join("".join(parse_subscript(s) for s in l) for l in args[1::2])

        if len(args) % 2 == 1:
            equation += "->" + "".join(parse_subscript(s) for s in args[-1])
            operands = args[:-1:2]
        else:
            operands = args[::2]
    else:
        equation = args[0]
        operands = args[1:]

    if has_torch_function(operands):
        return handle_torch_function(einsum, operands, equation, *operands)

    if len(operands) == 1 and isinstance(operands[0], (list, tuple)):
        _operands = operands[0]
        return einsum(equation, *_operands)

    if len(operands) <= 2 or not opt_einsum.enabled:
        return _VF.einsum(equation, operands)  # type: ignore[attr-defined]

    path = None
    if opt_einsum.is_available():
        key = (equation, tuple([op.shape for op in operands]))
        if key in path_dict:
            path = path_dict[key]
        else:
            _opt_einsum = opt_einsum.get_opt_einsum()
            tupled_path = _opt_einsum.contract_path(
                equation, *operands, optimize=opt_einsum.strategy
            )[0]
            path = [item for pair in tupled_path for item in pair]
            path_dict[key] = path
    return _VF.einsum(equation, operands, path=path)  # type: ignore[attr-defined]

def xop3(a1,a2,x):
    return a2 @ x @ a1
    
def xop5(a1,a2,b1,b2,x):
    return a2 @ x @ a1 + b2 @ x @ b1

xop3_vmap = pt.vmap(xop3,in_dims=(None,None,0))
xop5_vmap = pt.vmap(xop5,in_dims=(None,None,None,None,0))

def xop_fast(*args: Tensor) -> Tensor:
    if len(args) == 3:
        x = args[-1]
        
        y_shape = list(x.shape)
        y_shape[-1] = args[0].shape[1]
        y_shape[-2] = args[1].shape[0]
        
        return xop3_vmap(
            args[0],
            args[1],
            x.reshape((-1,x.shape[-2],x.shape[-1]))
        ).reshape(y_shape)
            
    elif len(args) == 5:
        x = args[-1]
        
        y_shape = list(x.shape)
        y_shape[-1] = args[0].shape[1]
        y_shape[-2] = args[1].shape[0]
        
        return xop5_vmap(
            args[0],
            args[1],
            args[2],
            args[3],
            x.reshape((-1,x.shape[-2],x.shape[-1]))
        ).reshape(y_shape)
            
    else:
        raise Exception(f'Error: expected 3 or 5 arguments to xop_fast, received {len(args)}')
        
def xop(
    M1: Tensor,
    M2: Tensor,
    x: Tensor
    ) -> Tensor:
#   tensor-product (m1 (x) m2) x

    shape = np.array(x.shape).tolist()
    n2 = shape[-2]
    n1 = shape[-1]
    
    y = x.reshape(-1,n2,n1) # add -3 index if it's missing
    
    if M1 is not None:
        M1_diag = False
        if len(M1.shape) == 2:
            m1,n1=M1.shape
            shape[-1] = m1
        else:
            M1_diag = True
        
    if M2 is not None:
        M2_diag = False
        if len(M2.shape) == 2:
            m2,n2=M2.shape
            shape[-2] = m2
        else:
            M2_diag = True
        
    if M1 is None:    
        if M2 is None:
            eq = 'ijk->ijk'
            ops = [y]
        else:
            if M2_diag:
                eq = 'j,ijk->ijk'
                ops = [M2, y]
            else:
                eq = 'lj,ijk->ilk'
                ops = [M2, y]
    else:
        if M2 is None:
            if M1_diag:
                eq = 'k,ijk->ijk'
                ops = [M1,y]
            else:
                eq = 'mk,ijk->ijm'
                ops = [M1,y]
        else:
            if M1_diag:
                if M2_diag:
                    eq = 'k,j,ijk->ijk'
                    ops=[M1,M2,y]
                else:
                    eq = 'k,nj,ijk->ink'
                    ops=[M1,M2,y]
            else:
                if M2_diag:
                    eq = 'mk,j,ijk->ijm'
                    ops=[M1,M2,y]
                else:
                    eq = 'mk,nj,ijk->inm'
                    ops=[M1,M2,y]
            
#   mx = pt.einsum(eq,*ops)
#   print(eq)
#   print(ops)
    mx = einsum(eq,*ops)
    return mx.reshape(shape)

class OpInv(pt.nn.Module):
    def __init__(self, Ax, Bx, Ay, By, dtype=None, device=None):
        super().__init__()
        
        if dtype is None and device is None:
            dtype  = Ax.dtype
            device = Ax.device
        
        args = { 'dtype': dtype, 'device': device }
        
        xdim = Ax.shape[0]
        ydim = Ay.shape[0]

        Ax = cf64(Ax)
        Bx = cf64(Bx)
        Ay = cf64(Ay)
        By = cf64(By)
        
        lamx,Sx = sp.linalg.eigh(Ax,Bx)
        lamy,Sy = sp.linalg.eigh(Ay,By)
        
        tol = 1e-4
        rank_def = False
        if np.abs(lamx[0]/lamx[1]) < tol and np.abs(lamy[0]/lamy[1]) < tol:
            rank_def = True
        
        self.Sx = Buffer(pt.tensor(Sx,**args))
        self.Sy = Buffer(pt.tensor(Sy,**args))
        
        Ix = pt.ones(xdim,dtype=pt.float64,device=pt.device('cpu'))
        Iy = pt.ones(ydim,dtype=pt.float64,device=pt.device('cpu'))
        
        lamx = pt.tensor(lamx,dtype=pt.float64,device=pt.device('cpu'))
        lamy = pt.tensor(lamy,dtype=pt.float64,device=pt.device('cpu'))
        
        diag = pt.kron(Iy,lamx) + pt.kron(lamy,Ix)
        Dinv=(1 / diag).reshape(ydim,xdim)
        
        if rank_def: # apply pseudo-inverse
            Dinv[0,0] = 0
        
        self.Dinv = Buffer(Dinv.to(**args))
        
    def forward(self, x: Tensor):
        shape = x.shape
        if len(shape) == 2:
            x.unsqueeze(0)
        
        res = xop(self.Sx,self.Sy,self.Dinv*xop(self.Sx.T,self.Sy.T,x))
        return res.reshape(shape)
    
class Op(pt.nn.Module):
    def __init__(
        self,
        Nx: int,
        Ny: int,
        dom: list,
        bc: str,
    ) -> None:
        super().__init__()
        
        self.Mbxn,self.Mbyn,self.Dbxn,self.Dbyn,self.Xn,self.Mbxm,self.Mbym,self.Dbxm,self.Dbym,self.Xm,self. Mbxp,self.Mbyp,self.Dbxp,self.Dbyp,self.Xp,self.Jxnm,self.Jynm,self.Jxpn,self.Jypn,self.Rx,self.Ry = setops((Nx,Ny),dom,bc,fullmass=True)
        
        self.D1x = Buffer(self.Jxpn.T @ self.Mbxn @ self.Dbxn)
        self.D1y = Buffer(self.Jypn.T @ self.Mbyn)
        
        self.D2x = Buffer(self.Jxpn.T @ self.Mbxn)
        self.D2y = Buffer(self.Jypn.T @ self.Mbyn @ self.Dbyn)
        
        self.Abx = Buffer(self.Dbxn.T @ self.Mbxn @ self.Dbxn)
        self.Aby = Buffer(self.Dbyn.T @ self.Mbyn @ self.Dbyn)
        
        # set up inverse operators
        Ax = self.Rx @ self.Dbxn.T @ self.Mbxn @ self.Dbxn @ self.Rx.T
        Ay = self.Ry @ self.Dbyn.T @ self.Mbyn @ self.Dbyn @ self.Ry.T
        Bx = self.Rx @ self.Mbxn @ self.Rx.T
        By = self.Ry @ self.Mbyn @ self.Ry.T
        
        self.Ainv = OpInv(Ax,Bx,Ay,By)
        
        self.Bxi=Buffer(pt.linalg.inv(self.Rx@self.Mbxn@self.Rx.T))
        self.Byi=Buffer(pt.linalg.inv(self.Ry@self.Mbyn@self.Ry.T))
        
        self.Binv = lambda x: xop(self.Bxi,self.Byi,x)

        Bx=self.Jxpn.T@self.Mbxn@self.Rx.T@self.Bxi@self.Rx@self.Mbxn@self.Jxpn
        By=self.Jypn.T@self.Mbyn@self.Ry.T@self.Byi@self.Ry@self.Mbyn@self.Jypn

        Ax=self.Jxpn.T@self.Mbxn@self.Dbxn@self.Rx.T@self.Bxi@self.Rx@self.Dbxn.T@self.Mbxn@self.Jxpn
        Ay=self.Jypn.T@self.Mbyn@self.Dbyn@self.Ry.T@self.Byi@self.Ry@self.Dbyn.T@self.Mbyn@self.Jypn
        
        self.Einv = OpInv(Ax,Bx,Ay,By)
        self.x = Buffer(self.Xn)
        self.last_refine = None
        
        self.JxnmT_Mbxm = Buffer(self.Jxnm.T@self.Mbxm)
        self.JynmT_Mbym = Buffer(self.Jynm.T@self.Mbym)
        
        self.Jxnm_Dbxn = Buffer(self.Jxnm @ self.Dbxn)
        self.Jynm_Dbyn = Buffer(self.Jynm @ self.Dbyn)
        
        self.RT_Binv_R_DT_1x = Buffer(self.Rx.T @ self.Bxi @ self.Rx @ self.D1x.T)
        self.RT_Binv_R_DT_2x = Buffer(self.Rx.T @ self.Bxi @ self.Rx @ self.D2x.T)
        self.RT_Binv_R_DT_1y = Buffer(self.Ry.T @ self.Byi @ self.Ry @ self.D1y.T)
        self.RT_Binv_R_DT_2y = Buffer(self.Ry.T @ self.Byi @ self.Ry @ self.D2y.T)
        
        self.RT_Binv_R_DT = lambda x: pt.stack((
            xop(self.RT_Binv_R_DT_1x,self.RT_Binv_R_DT_1y,x),
            xop(self.RT_Binv_R_DT_2x,self.RT_Binv_R_DT_2y,x)
        ))
        
        RT_S_x = Buffer(self.Rx.T @ self.Ainv.Sx)
        RT_S_y = Buffer(self.Ry.T @ self.Ainv.Sy)
        
        ST_R_x = Buffer(self.Ainv.Sx.T @ self.Rx)
        ST_R_y = Buffer(self.Ainv.Sy.T @ self.Ry)
        
        self.D1_RT_S_x = Buffer(self.D1x @ RT_S_x)
        self.D2_RT_S_x = Buffer(self.D2x @ RT_S_x)
        self.D1_RT_S_y = Buffer(self.D1y @ RT_S_y)
        self.D2_RT_S_y = Buffer(self.D2y @ RT_S_y)
        
        self.ST_R_D1T_x = Buffer(ST_R_x @ self.D1x.T)
        self.ST_R_D1T_y = Buffer(ST_R_y @ self.D1y.T)
        self.ST_R_D2T_x = Buffer(ST_R_x @ self.D2x.T)
        self.ST_R_D2T_y = Buffer(ST_R_y @ self.D2y.T)
        
        def D_RT_Ainv_R_DT_fn(x):
            t1 = pt.stack((
                xop(self.ST_R_D1T_x,self.ST_R_D1T_y,x),
                xop(self.ST_R_D2T_x,self.ST_R_D2T_y,x)
            ))
            
            t2 = self.Ainv.Dinv*t1
            return xop(self.D1_RT_S_x,self.D1_RT_S_y,t2[0]) + xop(self.D2_RT_S_x,self.D2_RT_S_y,t2[1])
        self.D_RT_Ainv_R_DT = D_RT_Ainv_R_DT_fn
        
    def Mb(self, x: Tensor):
        return xop_fast(self.Mbxn,self.Mbyn,x)
#       return xop(self.Mbxn,self.Mbyn,x)
        
    def Ab(self, x: Tensor):
        return xop_fast(self.Abx,self.Mbyn,self.Mbxn,self.Aby,x)
#       return xop(self.Abx,self.Mbyn,x) + xop(self.Mbxn,self.Aby,x)
        
    def Cb(self,cn,un):
        shape = un.shape
        
        vn = un.reshape(-1,shape[-2],shape[-1])
    
        dun = pt.stack((
#           xop(self.Dbxn,None,vn),
#           xop(None,self.Dbyn,vn)
            xop(self.Jxnm_Dbxn,self.Jynm,vn),
            xop(self.Jxnm,self.Jynm_Dbyn,vn)
        ))
        
        cu = xop(
            self.JxnmT_Mbxm,
            self.JynmT_Mbym,
            einsum(
                'ijk,injk->njk',
                xop(self.Jxnm,self.Jynm,cn),
                dun
            )
        )
        
        return cu.reshape(shape)
        
    def set_Hb(self, alpha=None) -> bool:
        if isinstance(alpha,list) and len(alpha) == 2:
            Abx = self.Dbxn.T @ self.Mbxn @ self.Dbxn
            Aby = self.Dbyn.T @ self.Mbyn @ self.Dbyn
            self.Hbx=Buffer((alpha[0]*.5)*self.Mbxn+alpha[1]*Abx)
            self.Hby=Buffer((alpha[0]*.5)*self.Mbyn+alpha[1]*Aby)
                
            Hx = self.Rx @ self.Hbx @ self.Rx.T
            Hy = self.Ry @ self.Hby @ self.Ry.T
            Bx = self.Rx @ self.Mbxn @ self.Rx.T
            By = self.Ry @ self.Mbyn @ self.Ry.T
            
            self.Hinv = OpInv(Hx,Bx,Hy,By)
            
            self.RxTSx_H = Buffer(self.Rx.T @ self.Hinv.Sx)
            self.RyTSy_H = Buffer(self.Ry.T @ self.Hinv.Sy)
            
            self.SxTRx_H = Buffer(self.Hinv.Sx.T @ self.Rx)
            self.SyTRy_H = Buffer(self.Hinv.Sy.T @ self.Ry)
            
            def RT_Hinv_R_fn(x):
                t1=xop(self.SxTRx_H,self.SyTRy_H,x)
                t2=self.Hinv.Dinv * t1
                t3 = xop(self.RxTSx_H,self.RyTSy_H,t2)
                return t3
                
            self.RT_Hinv_R = RT_Hinv_R_fn
            
            return True
        else:
            return False
        
    def Hb(self, x: Tensor, alpha=None):
        self.set_Hb(alpha)
#       return xop(self.Hbx,self.Mbyn,x) + xop(self.Mbxn,self.Hby,x)
        return xop_fast(self.Hbx,self.Mbyn,self.Mbxn,self.Hby,x)

    def Db(self, x: Tensor):
        return xop(self.D1x,self.D1y,x[0]) + xop(self.D2x,self.D2y,x[1])
        
    def R(self, x: Tensor):
        return xop(self.Rx,self.Ry,x)
        
    def RT(self, x: Tensor):
        return xop(self.Rx.T,self.Ry.T,x)
        
    def D(self, x: Tensor, restrict=False):
            y0, y1 = x[0], x[1]
            if restrict:
                y0 = self.RT(y0)
                y1 = self.RT(y1)
                
            return xop(self.D1x,self.D1y,y0) + xop(self.D2x,self.D2y,y1)
            
    def DT(self, x: Tensor, restrict=False):
        res = pt.stack((
            xop(self.D1x.T,self.D1y.T,x),
            xop(self.D2x.T,self.D2y.T,x)
        ))
        if restrict:
            res = self.R(res)
        return res
        
    def S(self, x: Tensor):
        xten = pt.tensor(x)
        xtype = xten.dtype
        nx = self.Jxpn.shape[1]
        ny = self.Jypn.shape[1]
        xten = xten.reshape(ny,nx).to(dtype=self.Jxpn.dtype)
        Sx = self.D_RT_Ainv_R_DT(xten)
#       Sx = self.D(self.Ainv(self.DT(xten,restrict=True)),restrict=True)
        o = Sx.flatten().to(dtype=xtype).numpy()
        return o
        
    def mask(self, x: Tensor):
        return self.RT(self.R(x))
        
    def maskc(self, x: Tensor):
        return x - self.mask(x)
        
    def incomp(self,usf,ps,bdti1):
        dp=-ortho(self.Einv(self.D(usf)));
        p=ps+dp*bdti1;
        uf=usf+self.RT_Binv_R_DT(dp)
    
        return uf,p
        
    def zero(self,var):
        if var == 'u':
            x_dim, y_dim = self.Mbxn.shape[0], self.Mbyn.shape[0]
            return pt.zeros(2, y_dim, x_dim)
                
        if var == 'p':
            x_dim, y_dim = self.Jxpn.shape[1], self.Jypn.shape[1]
        else: # scalar in mesh 1
            x_dim, y_dim = self.Mbxn.shape[0], self.Mbyn.shape[0]
            
        return pt.zeros(y_dim, x_dim)
        
    def refine(self,u,p,T=None,Nr=(1024,1024)):
        if self.last_refine is None or self.last_refine != Nr:
            xr0 = pt.linspace(self.x[0].min(),self.x[0].max(),Nr[0])
            yr0 = pt.linspace(self.x[1].min(),self.x[1].max(),Nr[1])
            
            self.Jxnl = J(xr0, self.x[0]).to(self.x[0].device,dtype=self.x[0].dtype)
            self.Jynl = J(yr0, self.x[1]).to(self.x[1].device,dtype=self.x[1].dtype)
            
        ur = xop(self.Jxnl,self.Jynl,u)
        pr = xop(self.Jxnl @ self.Jxpn, self.Jynl @ self.Jypn, p)
        if T is not None:
            tr = xop(self.Jxnl,self.Jynl,T)
        else:
            tr = None
       
        self.last_refine = Nr
        return ur, pr, tr
