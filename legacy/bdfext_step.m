function [un,bdti1]=bdfext_step(u,ub,hinv,hop,hop2,Bx,By,c,p,t,dt,istep)
    persistent ulag clag plag tlag;

    n=length(u);
    if (istep == 1) ulag=zeros(n,3); clag=ulag; plag=ulag; tlag=zeros(1,3); end

    ulag(:,3)=ulag(:,2);
    ulag(:,2)=ulag(:,1);
    ulag(:,1)=u;

    ic=min(3,istep);
    [betas,alphas]=bdfext(ic);
    disp(sum(betas));
    bdti=betas./dt;
    disp(alphas);
    b=-ulag(:,1)*bdti(2);
    if (istep>1); b=b-ulag(:,2)*bdti(3); end
    if (istep>2); b=b-ulag(:,3)*bdti(4); end
    b=a2u(By,Bx,b);

    disp(bdti(1))
    bdti1=bdti(1);
    un=(1./bdti1)*a2u(inv(By),inv(Bx),b);
end
