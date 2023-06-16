function [usf,ps,bdti1]=advance(uf,ub,p,t,dt,istep,Hb,Hinv,Bb,Binv,DbT,R,RT,C)
    persistent ulag clag plag tlag;

    c=C(uf,uf);

    if (istep == 1)
        n=length(uf); np=length(p); nc=length(c);
        ulag=zeros(n,3);
        clag=zeros(nc,3);
        plag=zeros(np,3);
        tlag=zeros(1,3);
    end

    ulag=lag(ulag,uf); clag=lag(clag,c); plag=lag(plag,p); tlag=lag(tlag,t);

    [bdti,alphas1,alphas2]=bdfext_var([t+dt tlag(1:min(3,istep))]',3);
    bdti1=bdti(1);

    b=-Bb(ulag*bdti(2:end));
    b=b-clag*alphas1;
    b=b-Hb(ub,bdti(1));
    ps=plag*alphas2;
    b=R(b+DbT(ps));

    usf=RT(Hinv(b,bdti1,istep<4))+ub;
end
