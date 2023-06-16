dt=5e-3;
T=100;
iostep=round(1/dt);
nsteps=round(T/dt);
iostep=nsteps;

conv_time=0.;
step_time=0.;
post_time=0.;

ce_time=0.;
ca_time=0.;

n=round(length(uf)*.5);
u=proj(uf,ub,Bb,nb);
ic=u;

blag=zeros(nb+1,3);
clag=zeros(nb+1,3);

vf=recon(u,ub);
t=0;
pm(vf); title(['Time=',num2str(t)]);

maxT=maxNumCompThreads(1);
c0=zeros(nb+1,1);

for istep=1:nsteps
    tic;
    k=min(istep,3);
    [betas,alphas]=bdfext(k);

%   c=reshape(reshape(Cbr,(nb+1)^2,nb+1)*u,nb+1,nb+1)*u;

    tic;

    c=c00+C_sub*u(2:end);

    tic;
    ce=reshape(reshape(Cjk,nb*(nb+1),nb)*u(2:end),nb+1,nb)*u(2:end);
    ce_time=ce_time+toc;

    tic;
    ca=Cp*((un2*u).*(uxn2*u)+(vn2*u).*(uyn2*u));
    ca_time=ca_time+toc;

    if ifexactc
        c=c+ce;
    else
        c=c+ca;
    end

    cerr=max(abs(ce-ca))/max(abs(ce));

    conv_time=conv_time+toc;

    b=Bbr*u;

    blag=lag(blag,b);
    clag=lag(clag,c);

    if istep<4;
        Hbr=nu*Abr+(betas(1)/dt)*Bbr; Hbr=.5*(Hbr+Hbr');
        Hbri=inv(Rr*Hbr*Rr');
    end

    r=-blag*(betas(2:end)'/dt);
    r=r-clag*alphas';
    r=r-Hbr(:,1);
    u(2:end)=Hbri*(Rr*r);

    t=t+dt;

    step_time=step_time+toc;

    if mod(istep,iostep) == 0
        tic;
        disp(['istep=',num2str(istep),' Cerr=',num2str(cerr)]);
        vf=recon(u,ub);
        pm(vf); title(['Time=',num2str(t)]);
        drawnow
        post_time=post_time+toc;
    end
end

vf=recon(u,ub);

disp(['conv_time=',num2str(conv_time)]);
disp(['step_time=',num2str(step_time)]);
disp(['post_time=',num2str(post_time)]);

disp(['ce_time=',num2str(ce_time)]);
disp(['ca_time=',num2str(ca_time)]);

maxNumCompThreads(maxT);
