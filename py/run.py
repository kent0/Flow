from op2d import Op
import cProfile
import io
from setops import setops
from problem import problem
from advance import BDFEXT, Solution
import matplotlib.pyplot as plt
import torch as pt
from util import mag
from vis import vis
import numpy as np
from ns import ns
from stokes import stokes
import pstats
import os

def main():
    Nx=128
    Ny=Nx
    tf=2e2
#   tf=4e1
    dt=2.5e-3
    nu=1/3000
    iotime=tf / 1000
    iptime=1
    
    nuis = pt.linspace(3000,4000,21)
    nuis = nuis[1:]
    
    pt.set_default_dtype(pt.float64)

    cname='ldc_reg_ref'
    pt.set_printoptions(sci_mode=True,precision=4)

    sol_stokes,op = stokes(Nx,Ny,cname,heat=True)
    vis(sol_stokes,op,folder='stokes/')
    
    for nui in nuis:
        nu = 1 / nui
        nui_str = str(int(pt.round(nui)))
        dir = f'data/{nui_str}/'
        os.makedirs(dir, exist_ok=True)
        sol,op = ns(Nx,Ny,tf,dt,cname,nu,u0=sol_stokes,iptime=iptime,iotime=iotime,pfx=dir)

if __name__ == "__main__":
#   cProfile.run('main()', sort='time')
    
    profile=False
    if profile:
        pr = cProfile.Profile()
        pr.run('main()')
        s = io.StringIO()
        ps = pstats.Stats(pr, stream=s)
        
        ps.sort_stats('time')
        ps.print_stats(20)
        
        ps.sort_stats('cumulative')
        ps.print_stats(20)
        
        print(s.getvalue())
    else:
        main()