function [un,ps,bdti1]=bdfext_step2(u,ub,hinv,hop,mass,DbT,c,p,t,dt,istep)
    persistent ulag clag plag tlag;

    if (istep == 1)
        n=length(u); np=length(p);
        ulag=zeros(n,3);
        clag=ulag;
        plag=zeros(np,3);
        tlag=zeros(1,3);
    end

    ulag=lag(ulag,u); clag=lag(clag,c); plag=lag(plag,p); tlag=lag(tlag,t);

    [bdti,alphas1,alphas2]=bdfext_var([t+dt tlag(1:min(3,istep))]',3);
    b=-mass(ulag*bdti(2:end));

    b=b-clag*alphas1;
    ps=plag*alphas2;
    b=b+DbT(ps);

    b=b-hop(ub,bdti(1));

    bdti1=bdti(1);

    un=hinv(b,bdti1,istep<4);
end
