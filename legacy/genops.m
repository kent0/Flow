Abr=zeros(nb+1,nb+1);
Bbr=zeros(nb+1,nb+1);
Cbr=zeros(nb+1,nb+1,nb+1);
C_approx=zeros(nb+1,nb+1);

Rr=speye(nb+1); Rr=Rr(2:end,:);

for j=1:nb+1
    Auj=Ab(ub(:,j));
    Buj=Bb(ub(:,j));
    for i=1:nb+1
        Abr(i,j)=dot(ub(:,i),Auj);
        Bbr(i,j)=dot(ub(:,i),Buj);
    end
end

Abr=.5*(Abr+Abr');
Bbr=.5*(Bbr+Bbr');

for k=1:nb+1
    disp(['Gen C k=',num2str(k)]);
    for j=1:nb+1
        Cujk=C(ub(:,k),ub(:,j));
        for i=1:nb+1
            Cbr(i,j,k)=dot(ub(:,i),Cujk);
        end
    end
end

c00=Cbr(:,1,1);
C0k=Cbr(:,2:end,1); C0k=reshape(C0k,nb+1,nb);
Cj0=Cbr(:,1,2:end); Cj0=reshape(Cj0,nb+1,nb);
Cjk=Cbr(:,2:end,2:end);

Juv=a2u(Jynm,Jxnm,ub);

BJx=Bbxm*Jxn2m;
BJy=Bbym*Jyn2m;

cbm=a2u(Jyn2m,Jxn2m,cb_vec);

% Deim
nc=nb+1;
%nc=nb*2;

for j=1:nc
%   Bcj=a2u(BJy,BJx,cb_vec(:,j));
    Bcj=a2u(Bbym,Bbxm,cbm(:,j));
    for i=1:nb+1
        C_approx(i,j)=dot(Juv(:,i),Bcj);
    end
end

[PUi,ps,U]=deim(cb_vec(1:n2,:),nc);

un2=a2u(Jynn2,Jxnn2,ub);
uxn2=a2u(Jynn2,Jxnn2,a2u(Iby,Dbx,ub));
uyn2=a2u(Jynn2,Jxnn2,a2u(Dby,Ibx,ub));

vn2=un2(ps'+n2,1:nb+1)*Rr'*Rr;
un2=un2(ps',1:nb+1)*Rr'*Rr;

uxn2=uxn2(ps',1:nb+1)*Rr'*Rr;
uyn2=uyn2(ps',1:nb+1)*Rr'*Rr;

C_sub=C0k+Cj0;
Cp=C_approx*PUi;
