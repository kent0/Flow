function [ee,s,S,it,pmat,bmat,amat,c1mat,c2mat]=sadvance2(uf,S,Hb,Hinv,Bb,Binv,Einv,Db,DT,R,RT,C,pm,ptype)

    tol=1.e-8;
%   tol=1.e-15;
    restart=25;
    maxit=2000;
    eps=0;

    [b1,b2]=evalf(uf,C,Hb,Binv,Einv,Db,DT,R,RT);
    b=b1+b2;
    ee=norm(b);

    x=b*0;
    aop=@(s) R(Hb(RT(s),0));
    c1op=@(s) R(C(RT(s),uf));
    c2op=@(s) R(C(uf,RT(s)));
    bop=Binv;
    acop=@(s) Binv(R(Hb(RT(s),0)+C(RT(s),uf)+C(uf,RT(s))));
    pop=@(s)R(incomp2(RT(s),Binv,Einv,Db,DT,RT));
    jop=@(s) pop(acop(s));

    amat=[];
    bmat=[];
    c1mat=[];
    c2mat=[];
    pmat=[];

%   acmat=opmat(acop,length(x));
%   jmat=opmat(jop,length(x));
%   amat=opmat(aop,length(x));
%   c1mat=opmat(c1op,length(x));
%   c2mat=opmat(c2op,length(x));
%   bmat=opmat(bop,length(x));
%   pmat=opmat(pop,length(x));

%   e1=eig(acmat);
%   e2=eig(jmat);
%   e3=eig(bmat*amat);

%   figure; scatter(real(e1),imag(e1)); hold on; scatter(real(e2),imag(e2),'*');
%   scatter(real(e3),imag(e3),'+');
%   ellipse(e1);
%   ellipse([e2;max(abs(e3))]);
%   size(e2)

    b=pop(b);
    bnorm=norm(b)/norm(uf)
        btol=100;

    if ptype == 1
        minv=@(s)s;
    elseif ptype == 2
        minv=@(s) pcr2(s,acop,iden,maxit,1.e-10,[]);
    elseif ptype == 3
        if bnorm > btol
            minv=@(s) 0.0005*polyinv(@(s)acop(s*0.0005),s,20);
        else
            minv=@(s) s;
        end
    elseif ptype == 4
%       acmat=opmat(acop,length(x));
        amat=opmat(aop,length(x));
        c1mat=opmat(c1op,length(x));
        c2mat=opmat(c2op,length(x));
        ac=amat+c1mat+c2mat;
        im=inv(ac);
        minv=@(s) im*s;
    elseif ptype == 5
        m=50;
%       scal=1+20*exp(-bnorm);
%       scal=1/exp(-bnorm/20)
%       scal=1;

        scal=20;
        mop=@(s) R(scal*Hb(RT(s),0)+C(uf,RT(s)));
        mop=@(s) R(Hb(RT(s),0));
        [Q,H]=arnoldi(mop,b,m);

        es=eig(H(1:m,1:m));
        [a,bb,c,x0]=ellipse(es);

        d=x0;
        c=c*1;

        niter=20;
        rat=exp(-(bnorm-btol)/10000)
        if bnorm > btol
            minv=@(s) (1-rat)*cheby3(mop,s,c,d,niter)+rat*s;
        else
            minv=@(s) s;
        end
    end

    s0=b*0;
    if length(S) > 0
        as=S*0;
        for i=1:size(S,2)
            as(:,i)=jop(S(:,i));
        end
        gamma=(as'*as)\(as'*-b);
        s0=S*gamma;
    end
%   s0=s0*0;

%   m=20;
%   [Q,H] = arnoldi(jmat,-b,m);
%   es=eig(H(1:m,1:m));
%   hold on; scatter(real(es),imag(es),'^')

    [s,e,it,es]=mygmres(jop,minv,-b,s0,maxit,tol);
%   hold on; scatter(real(es),imag(es),'x')
%   [s,e,it,es]=mygmres(jop,minv,rand(length(b),1),s0,maxit,tol);
%   hold on; scatter(real(es),imag(es),'v')
%   pause
    s=pop(s);
    S=orthoadd(S,s);
    it
end
