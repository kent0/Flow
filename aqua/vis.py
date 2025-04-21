from matplotlib import pyplot as plt
from matplotlib.image import imsave
from matplotlib.cm import get_cmap
import torch as pt
from op2d import Op
from util import mag

def vis(sol,op: Op,T=None,folder='cache/',i=None):
    
    if i is None:
        istr = ''
    else:
        istr = f'{i}'.rjust(6,'0')
        
    u = sol.u
    p = sol.p
    T = sol.T
    
    v,q,S = op.refine(u,p,T)
    
    v = v.flip(dims=(1,))
    vmag = mag(v)
    
    q = q.flip(dims=(0,))
    if S is not None:
        S = S.flip(dims=(0,))
    
    vmag = vmag.cpu()
    v = v.cpu()
    q = q.cpu()
    if S is not None:
        S = S.cpu()
    
    print(f'umag: min,max = {vmag.min().item():.3f},{vmag.max().item():.3f}')
    print(f'u: min,max = {v[0].min().item():.3f},{v[0].max().item():.3f}')
    print(f'v: min,max = {v[1].min().item():.3f},{v[1].max().item():.3f}')
    print(f'p: min,max = {q.min().item():.3f},{q.max().item():.3f}')
    if S is not None:
        print(f'T: min,max = {S.min().item():.3f},{S.max().item():.3f}')
    
    cm = get_cmap('turbo',2048)
    
    imsave(folder+f'umag{istr}.png', vmag,cmap=cm)
    imsave(folder+f'u{istr}.png', v[0],cmap=cm)
    imsave(folder+f'v{istr}.png', v[1],cmap=cm)
    imsave(folder+f'p{istr}.png', q,cmap=cm)
    if S is not None:
        imsave(folder+f'T{istr}.png', S,cmap=cm)
        