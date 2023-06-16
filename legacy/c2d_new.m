function cu=c2d_new(cxyn,uvn,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jx,Jy,Fx,Fy,Qx,Qy,Rx,Ry)
    persistent um vm uxm uym vxm vym cxm cym t cu;

    udim=2;
    nn=round(length(uvn)/udim);
    um=a2u(Jy,Jx,uvn(1:nn));
    vm=a2u(Jy,Jx,uvn(nn+1:end));

    uxm=a2u(Iby,Dbx,um);
    uym=a2u(Dby,Ibx,um);
    vxm=a2u(Iby,Dbx,vm);
    vym=a2u(Dby,Ibx,vm);

    cxm=a2u(Jy,Jx,cxyn(1:nn));
    cym=a2u(Jy,Jx,cxyn(nn+1:end));

    t=[cxm.*uxm+cym.*uym; ...
       cxm.*vxm+cym.*vym];

    cu=a2u(Ry,Rx,a2u(Qy',Qx',a2u(Jy',Jx',a2u(Bby,Bbx,t))));
end
