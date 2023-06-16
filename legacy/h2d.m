function hu=h2d(ub,Abx,Aby,Bbx,Bby,nu,bdti)
    Hx=nu*Abx+(bdti*.5)*Bbx;
    Hy=nu*Aby+(bdti*.5)*Bby;
    hu=a2u(Bby,Hx,ub)+a2u(Hy,Bbx,ub);
end
