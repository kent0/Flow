function Juv=interp_uv(uv,Nx,Ny,Mx,My);
    zxn=zwgll(Nx);
    zyn=zwgll(Ny);
    zxm=zwgll(Mx);
    zym=zwgll(My);
    Jx=interp_mat(zxm,zxn);
    Jy=interp_mat(zym,zyn);
    Juv=a2u(Jy,Jx,uv);
end
