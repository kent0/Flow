A_rom=zeros(nb+1,nb+1);
B_rom=A_rom;
C_rom=zeros(nb+1,nb+1,nb+1);
C_approx=zeros(nb+1,nb+1);

z=zb;

for j=1:nb+1
    Azj=Ab(z(:,j));
    Bzj=Bb(z(:,j));
    for i=1:nb+1
        Abr(i,j)=dot(z(:,i),Azj);
        Bbr(i,j)=dot(z(:,i),Bzj);
    end
end

for k=1:nb+1
    disp(['C gen k=',num2str(k)]);
    for j=1:nb+1
        Czjk=C(z(:,k),z(:,j));
        for i=1:nb+1
            Cbr(i,j,k)=dot(z(:,i),Czjk);
        end
    end
end

I=speye(nb+1);
Rr=I(2:nb+1,:);

deim;

Jz=a2u(Jynm,Jxnm,z);

Cij0=zeros(nb+1,nb+1);
Ci0k=zeros(nb+1,nb+1);
Ci00=zeros(nb+1,1);

C00=C(z(:,1),z(:,1));

for j=1:nb+1
    Cj0=C(z(:,1),z(:,j));
    C0k=C(z(:,j),z(:,1));
    Ci00(j)=dot(C00,z(:,j));
    for i=1:nb+1
        Cij0(i,j)=dot(Cj0,z(:,i));
        Ci0k(i,j)=dot(C0k,z(:,i));
    end
end

for j=1:nb+1
    Bcj=a2u(Bbym,Bbxm,csnap(:,j));
    for i=1:nb+1
        C_approx(i,j)=dot(Jz(:,i),Bcj);
    end
end
