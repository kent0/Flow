function cfl=gcfl(uf,dt,Xn,Yn);

cflx=0.;
cfly=0.;

n2=round(length(uf)/2);
[Nx,Ny]=size(Xn); Nx=Nx-1; Ny=Ny-1;
Uf=reshape(uf(1:n2),Nx+1,Ny+1);
Vf=reshape(uf(n2+1:end),Nx+1,Ny+1);

for j=1:Ny+1
for i=1:Nx+1
    if (i~=1 && i~=Nx+1)
        dxmin=min(Xn(i+1,j)-Xn(i,j),Xn(i,j)-Xn(i-1,j));
    elseif (i==1)
        dxmin=Xn(i+1,j)-Xn(i,j);
    elseif (i==Nx+1)
        dxmin=Xn(i,j)-Xn(i-1,j);
    end
    if (j~=1 && j~=Nx+1)
        dxmin=min(Xn(i,j+1)-Xn(i,j),Xn(i,j)-Xn(i,j-1));
    elseif (i==1)
        dxmin=Xn(i,j+1)-Xn(i,j);
    elseif (i==Nx+1)
        dxmin=Xn(i,j)-Xn(i,j-1);
    end
    cflxi=Uf(i,j)*dt/dxmin; cflx=max(cflxi,cflx);
    cflyi=Vf(i,j)*dt/dymin; cfly=max(cflyi,cfly);
end
end

cfl=max(cflx,cfly);
