function [ee,S]=sadvance(uf,ub,S,Hb,Hinv,Bb,Binv,Einv,Db,DT,R,RT,C,sigma,pm)

    tol=1.e-9;
    restart=25;
    maxit=2000;

    [b1,b2]=evalf(uf,C,Hb,Binv,Einv,Db,DT,R,RT);
    b=b1+b2;
    ee=norm(b);
    eps=1e-6;
    pm(RT(b1))
    figure
    pm(RT(b2))
    figure
    pm(incomp2(RT(b2),Binv,Einv,Db,DT,RT))
    pause

%   b=R(incomp2(RT(Binv(R(Hb(uf,0)+C(uf,uf)))),Binv,Einv,Db,DT,RT));
%   e=norm(b)/norm(uf)

    x=b*0;
%   if length(s) > 2; x=s; end
    jop=@(s) R(incomp2(RT(Binv(R(Hb(RT(s),0)+C(RT(s),uf)+C(uf,RT(s))))),Binv,Einv,Db,DT,RT)) + eps*s;

    aop=@(s) R(Hb(RT(s),0));
    c1op=@(s) R(C(RT(s),uf));
    c2op=@(s) R(C(uf,RT(s)));
    bop=Binv;
    pop=@(s)R(incomp2(RT(s),Binv,Einv,Db,DT,RT));

    jmat=opmat(jop,length(x));
    amat=opmat(aop,length(x));
    c1mat=opmat(c1op,length(x));
    c2mat=opmat(c2op,length(x));
    bmat=opmat(bop,length(x));
    pmat=opmat(pop,length(x));

%   ej=eig(jmat);
%   figure
%   scatter(real(ej),imag(ej));
%   title('jmat');

%   e=eig(bmat*(amat+c1mat+c2mat));
%   figure
%   scatter(real(e),imag(e));
%   title('acmat');

%   ea=eig(amat);
%   figure
%   scatter(real(ea),imag(ea));
%   title('amat');

%   ec1=eig(c1mat);
%   figure
%   scatter(real(ec1),imag(ec1));
%   title('c1mat');

%   ec2=eig(c2mat);
%   figure
%   scatter(real(ec2),imag(ec2));
%   title('c2mat');

%   eb=eig(bmat);
%   figure
%   scatter(real(eb),imag(eb));
%   title('bmat');

%   ep=eig(pmat);
%   figure
%   scatter(real(ep),imag(ep));
%   title('pmat');

%   e1=eig((bmat*(amat))\jmat);
%   figure
%   scatter(real(e1),imag(e1));
%   title('e1');

%   e2=eig((bmat*(amat+c1mat))\jmat);
%   figure
%   scatter(real(e2),imag(e2));
%   title('e2');

%   e3=eig((bmat*(amat+c2mat))\jmat);
%   figure
%   scatter(real(e3),imag(e3));
%   title('e3');

%   e4=eig((bmat*(amat+c1mat+c2mat))\jmat);
%   figure
%   scatter(real(e4),imag(e4));
%   title('e4');

%   e5=eig(jmat*inv(bmat*(amat+c1mat+c2mat)));
%   figure
%   scatter(real(e5),imag(e5));
%   title('e5');

    jjmat=(bmat*(amat+c1mat+c2mat))\jmat;
    jjop=@(s)jjmat*s;
    bb=(bmat*(amat+c1mat+c2mat))\b;
    minvmat=inv(bmat*(amat+c1mat+c2mat));

    iden=@(x)x;
    minv=iden;
    minv=@(x)minvmat*x;

%   figure
%   scatter(real(e-ej),imag(e-ej));

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

    [s,e,it]=mygmres(jop,minv,-b,s0,maxit,tol);
    s=R(incomp2(RT(s),Binv,Einv,Db,DT,RT));
    S=[S s];
    it
    usf=uf+sigma*RT(s);
end
