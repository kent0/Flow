function du=d2d(uv,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,ifrt)
    if (ifrt) rtuv=a2u(Ry',Rx',uv); else rtuv=uv; end
    n2=round(length(rtuv)/2);
    u=rtuv(1:n2); v=rtuv(n2+1:end);
%   du=a2u(Jpy',Jpx',a2u(Bby,Bbx,(a2u(Iby,Dbx,u)+a2u(Dby,Ibx,v))));
    du=a2u(Jpy'*Bby,Jpx'*Bbx*Dbx,u)+a2u(Jpy'*Bby*Dby,Jpx'*Bbx,v);
end

%function du=d2d(uv,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,ifrt)
%    if (ifrt) rtuv=a2u(Ry',Rx',uv); else rtuv=uv; end
%    n2=round(length(rtuv)/2);
%    u=rtuv(1:n2); v=rtuv(n2+1:end);
%    du=a2u(Jpy',Jpx',a2u(Bby,Bbx,(a2u(Iby,Dbx,u)+a2u(Dby,Ibx,v))));
%end
