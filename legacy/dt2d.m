function dtp=dt2d(p,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,ifr)
    t=a2u(Bby,Bbx,a2u(Jpy,Jpx,p));
%   t=a2u(Jpy,Jpx,p);
    dtp=[a2u(Iby,Dbx',t);a2u(Dby',Ibx,t)];
    if ifr; dtp=a2u(Ry,Rx,dtp); end
end
