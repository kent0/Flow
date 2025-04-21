import torch as pt
from torch import Tensor
from torch import nn
from torch.nn import Buffer
from typing import Optional

class Solution(nn.Module):
    u: Optional[Tensor]
    p: Optional[Tensor]
    T: Optional[Tensor]
    def __init__(self, u, p=None, T=None):
        super().__init__()
        self.u = Buffer(u) if u is not None else None
        self.p = Buffer(p) if p is not None else None
        self.T = Buffer(T) if T is not None else None
        

def cf64(x: Tensor):
    return x.cpu().to(dtype=pt.float64)

def mag(u):
    return pt.sqrt(u[0]**2 + u[1]**2)

