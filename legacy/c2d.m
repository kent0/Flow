function cu=c2d(cn,uvn,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jx,Jy,Fx,Fy,Rx,Ry)
    udim=2;
    uvm=a2u(Jy,Jx,uvn);
    cm=a2u(Jy,Jx,cn);
    n=round(length(uvm)/udim);

    ux=a2u(Iby,Dbx,uvm);
    uy=a2u(Dby,Ibx,uvm);

    if udim == 2
        t=[cm(1:n).*ux(1:n)+cm(1+n:end).*uy(1:n); ...
           cm(1:n).*ux(1+n:end)+cm(1+n:end).*uy(1+n:end)];
    else
        t=[cm(1:n).*ux(1:n)+cm(1+n:end).*uy(1:n)];
    end

    cu=a2u(Jy',Jx',a2u(Bby,Bbx,t));
end
