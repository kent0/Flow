using Plots

function stokes(Nx,Ny,cname,pfx="")
    # usage: ns(32,32,100,1e-3,"ldc",1/1000,[]);

    ifconv=true;

    bc="dddd";

    nf=Int(round(Nx*.1));
    nf=0;
    df=.05;
    
    nu=0.5;
    nu=1;

    ifplot=true;
    #ifplot=false;

    (exact,f,dom,ifexact)=problem(cname,nu);

    Abx,Aby,Bbx,Bby,Bbxm,Bbym,Bbxn2,Bbyn2,Bbpx,Bbpy,Dbx,Dby,Dbxm,Dbym,Fx,Fy,Ibx,Iby,Ibxm,Ibym,Jxnm,Jynm,Jxmn,Jymn,Jxnn2,Jynn2,Jxn2m,Jyn2m,Jxn2l,Jyn2l,Jxml,Jyml,Jxnl,Jynl,Jxpn,Jypn,Jxpl,Jypl,Rx,Ry,Xn,Yn,Xn2,Yn2,Xm,Ym,Xl,Yl,Xpn,Ypn=setops(Nx,Ny,dom,bc,nf,df);

    x=reshape(Xn,:); y=reshape(Yn,:);
    xp=reshape(Xpn,:); yp=reshape(Ypn,:);

    (Ax,Ay,Bx,By,Dx,Dy,Ix,Iy)=rops(Abx,Aby,Bbx,Bby,Dbx,Dby,Rx,Ry);

    (Bxi,Byi,Bbpxi,Bbpyi)=invops(Bx,By,Bbpx,Bbpy);

    (Ab,B,Bb,Binv,Hb,Hinv,DbT,DT,Db,D,Einv,C,R,RT,
       vort,sf,maskc)=setfuns(Abx,Aby,Ax,Ay,
        Bbx,Bby,Bx,By,Bxi,Byi,Bbxm,Bbym,Dbx,Dby,Dbxm,Dbym,
        Ibx,Iby,Ix,Iy,Ibxm,Ibym,Jxnm,Jynm,Jxpn,Jypn,Fx,Fy,Rx,Ry,nu);

    if ~ifconv; C=(vf,uf) -> R(uf)*0; end

    pdim=2;

    (p1,p2,pm,pv,psf,pp,pdiv,plot_m1,plot_m2,plot_mm,plot_m3)=setplots(Xn,Yn,Xn2,Yn2,Xm,Ym,Xpn,Ypn,Jxnl,Jynl,Jxn2l,Jyn2l,Jxml,Jyml,Jxpl,Jypl,Xl,Yl,pdim,vort,sf,Db);

    global t=0;
    ue=exact(x,y,t);
    if cname == "kov"; pe=kov_p(xp,yp,nu); end
    if cname == "walsh"; pe=walsh_p(xp,yp,0,nu); end
    ub=maskc(ue);
    p=xp*0; dp=p;

    uf=ue;
    uf_old=uf;

    ttot=0.;
    count=0;

    ubsf=uf*0;

    ulag=[];
    clag=[];
    plag=[];
    tlag=[];

    dt = 1.e-6;
    ue_next=exact(x,y,t+dt); ub=maskc(ue_next);
    
    if bc=="pppp"; ub=ub*0; end

    Ab=(u,bdti)->Hb(u,0)
    Ainv=(b,bdti,ifupd)->Hinv(b,0,true)
    Ainv1=(b)->Hinv(b,0,true)
    
    Eainv=(p)->sa2dinv(p,Bbx,Bby,Dbx,Dby,Jxpn,Jypn,Rx,Ry);
    C_dummy=(cf,uf)->C(cf,uf)*0;
    
    istep=1
    
    uf=uf*0;
    usf,ps,bdti1=advance(uf,ub,p,t,dt,istep,Ab,Ainv,Bb,Binv,DbT,R,RT,C_dummy);
    bdti1=1;
    # should be solution to poisson
    pm(usf,"test")
    (uf,p)=incomp(usf,ps,bdti1,Ainv1,Eainv,Db,DT,RT);
    iostep=1
    nsteps=1

    (diff,err)=post(uf,ue,ue,p,ubsf,t,dt,Xn,Yn,istep,iostep,nsteps,
        ifexact,ifplot,Bb,p1,p2,pm,pv,psf,pp,pdiv,pfx);

    u=reshape(uf,Nx+1,Ny+1,2);

    println("Done!");
    return (u,err,ttot,dt)
end
