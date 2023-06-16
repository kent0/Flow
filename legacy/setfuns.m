function [Ab,B,Bb,Binv,Hb,Hinv,DbT,DT,Db,D,Einv,C,R,RT, ...
   vort,sf,maskc]=setfuns(Abx,Aby,Ax,Ay, ...
   Bbx,Bby,Bx,By,Bxi,Byi,Bbxm,Bbym,Dbx,Dby,Dbxm,Dbym, ...
   Ibx,Iby,Ix,Iy,Ibxm,Ibym,Jxm,Jym,Jpx,Jpy,Fx,Fy,Rx,Ry,nu);

Hinv=@(b,bdti,ifupd) h2dinv(b,Ax,Ay,Bx,By,Ix,Iy,nu,bdti,ifupd);
Hb=@(ub,bdti) h2d(ub,Abx,Aby,Bbx,Bby,nu,bdti);
Ab=@(uf) a2u(Bby,Abx,uf)+a2u(Aby,Bbx,uf);
B=@(u) a2u(By,Bx,u);
Bb=@(uf) a2u(Bby,Bbx,uf);
Binv=@(u) a2u(Byi,Bxi,u);

DT=@(p) dt2d(p,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,true);
DbT=@(p) dt2d(p,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,false);
D=@(uv) d2d(uv,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,true);
Db=@(uv) d2d(uv,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,false);

Einv=@(p) sb2dinv(p,Bbx,Bby,Dbx,Dby,Jpx,Jpy,Rx,Ry);

C=@(cf,uf) c2d(cf,uf,Bbxm,Bbym,Dbxm,Dbym,Ibxm,Ibym,Jxm,Jym,Fx,Fy,Rx,Ry);

R=@(uf) a2u(Ry,Rx,uf);
RT=@(u) a2u(Ry',Rx',u);

vort=@(uf) v2d(uf,Dbx,Dby);

RBbx=Rx*Bbx;
RBby=Ry*Bby;

RAbx=Rx*Abx;
RAby=Ry*Aby;

sf=@(b,ubsf) ubsf+a2u(Ry',Rx',a2dinv(a2u(RBby,RBbx,vort(b))-a2d(ubsf,RAby,RAbx,RBby,RBbx),Ax,Ay,Bx,By,Rx,Ry));

mask=@(u) a2u(Ry'*Ry,Rx'*Rx,u);
maskc=@(u) u-mask(u);
