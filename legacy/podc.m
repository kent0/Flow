n2=(2*Nx+1)*(2*Ny+1);
cs=zeros(n2*2,ndump);

GC_vec=zeros(ndump,ndump);
GC1=zeros(ndump,ndump);
GC2=zeros(ndump,ndump);

Jus=a2u(Jynn2,Jxnn2,us0);
Jusx=a2u(Jynn2,Jxnn2,a2u(Iby,Dbx,us0));
Jusy=a2u(Jynn2,Jxnn2,a2u(Dby,Ibx,us0));

for i=1:ndump
    cs(:,i)=[Jus(1:n2,i).*Jusx(1:n2,i)+Jus(n2+1:end,i).*Jusy(1:n2,i);
            Jus(1:n2,i).*Jusx(n2+1:end,i)+Jus(n2+1:end,i).*Jusy(n2+1:end,i)];
end

for j=1:ndump
    disp(['gengram C ',num2str(j)]);
    Bcj=a2u(Bbyn2,Bbxn2,cs(:,j));
    for i=1:ndump
        GC_vec(i,j)=dot(cs(:,i),Bcj);
        GC1(i,j)=dot(cs(1:n2,i),Bcj(1:n2));
        GC2(i,j)=dot(cs(n2+1:end,i),Bcj(n2+1:end));
    end
end

GC_vec=.5*(GC_vec+GC_vec');
GC1=.5*(GC1+GC1');
GC2=.5*(GC2+GC2');

[VC_vec,DC_vec]=eig(GC_vec);
[VC1,DC1]=eig(GC1);
[VC2,DC2]=eig(GC2);

DC_vec=flip(diag(DC_vec));
VC_vec=flip(VC_vec,2);

DC1=flip(diag(DC1));
VC1=flip(VC1,2);

DC2=flip(diag(DC2));
VC2=flip(VC2,2);

cb_vec=cs*VC_vec;

cb=[cs(1:n2,:)*VC1;cs(n2+1:end,:)*VC2];

for i=1:ndump
    cb_vec(:,i)=cb_vec(:,i)./sqrt(dot(cb_vec(:,i),a2u(Bbyn2,Bbxn2,cb_vec(:,i))));
    cb(1:n2,i)=cb(1:n2,i)./sqrt(dot(cb(1:n2,i),a2u(Bbyn2,Bbxn2,cb(1:n2,i))));
    cb(n2+1:end,i)=cb(n2+1:end,i)./sqrt(dot(cb(n2+1:end,i),a2u(Bbyn2,Bbxn2,cb(n2+1:end,i))));
end
