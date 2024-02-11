using Profile, PProf
include("flow.jl")

Nx=105;Ny=Nx;
T=100.0;
#T=0.1;
dt=3.125e-3;
cname="ldc_reg";
nu=1/(3e4);
#u0=[];
#u=ns(Nx,Ny,T,dt,cname,nu,u0);
@profile include("ns.jl")

pprof(;webport=58699)
