import torch as pt
import numpy as np

def ldc_reg_ref(x,y,t,nu):
#   yc = y.cpu()
#   xc = x.cpu()
    yc = y
    xc = x

    uc = np.stack((np.heaviside(y-(1-1e-3),0.5)*((1.0+x)*(1.0-x))**2,y*0),axis=0)
    return uc
    
def problem(cname,dtype=None,device=None):
    lam=pt.pi*0.5
    z = pt.tensor([0])
    pd = {
        'st2': {
            'ufun': lambda x,y,t,nu: st2_uv(x, y, t, w, nu),
            'dom': [0,1,0,1],
            'f': None,
            'exact': True
        },
        'start': {
            'ufun': lambda x,y,t,nu: start_uv(x,y,t,10),
            'dom': [-1,1,-1,1],
            'f': lambda x,y,t,nu: pt.stack((x*0+1.5*nu,y*0)),
            'exact': False
        },
        "kov":{
            'ufun': lambda x,y,t,nu: kov_uv(x,y,nu),
            'dom': [-0.5,2.0,-0.5,1.5],
            'f': None,
            'exact': True
        },
        "walsh": {
            'ufun': lambda x,y,t,nu: walsh_uv(x,y,t,nu),
            'dom': [0,2*pt.pi,0,2*pt.pi],
            'f': None,
            'exact': True
        },
        "ldc": {
            'ufun': lambda x,y,t,nu: pt.stack((pt.heaviside(y-(1-1e-3),z),y*0)),
            'dom': [0,1,0,1],
            'f': None,
            'exact': False
        },
        "ldc_reg": {
            'ufun': lambda x,y,t,nu: pt.stack((16*pt.heaviside(y-(1-1e-3),z)*(x*(1.0-x))**2,y*0)),
            'dom': [0,1,0,1],
            'f': None,
            'exact': False,
        },
        "ldc_reg_ref": {
            'ufun': ldc_reg_ref,
            'dom': [-1,1,-1,1],
            'f': None,
            'exact': False
        },
        "swirl": {
            'ufun': lambda x,y,t,nu: pt.stack((2*y*(1-x^2),-2*x*(1-y^2))),
            'dom': [-1,1,-1,1],
            'f': lambda x,y,t,nu: pt.stack((4*y,-4*x))*nu,
            'exact': False,
        },
        "decay": {
            'ufun': lambda x,y,t,nu: pt.stack((pt.cos(lam*y)*pt.exp(-lam^2*nu*t),x*0)),
            'dom': [-1,1,-1,1],
            'f': None,
            'exact': True
        },
        "decay2": {
            'ufun': lambda x,y,t,nu: pt.stack((1.5*(1-y*y),x*0))-start_uv(x,y,t,30),
            'dom': [-1,1,-1,1],
            'f': None,
            'exact': True
        }
    }
    return pd[cname]