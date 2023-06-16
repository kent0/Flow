function cfl=gcfl(uf,dt,Xn,Yn);
        cflx=0.;
        cfly=0.;
        n2=round(length(uf)/2);
        [Nx,Ny]=size(Xn); Nx=Nx-1; Ny=Ny-1;
        Uf=reshape(uf(1:n2),Nx+1,Ny+1);
        Vf=reshape(uf(n2+1:end),Nx+1,Ny+1);
        for j=2:Ny
        for i=2:Nx
            cflxi=Uf(i,j)*dt/min(Xn(i+1,j)-Xn(i,j),Xn(i,j)-Xn(i-1,j));
            if cflxi>cflx; cflx=cflxi; end
            cflyi=Vf(i,j)*dt/min(Yn(i,j+1)-Yn(i,j),Yn(i,j)-Yn(i,j-1));
            if cflyi>cfly; cfly=cflyi; end
        end
        end
        cfl=max(cflx,cfly);
