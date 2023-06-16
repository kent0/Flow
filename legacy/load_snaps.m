function us=load_snaps(Nx,Ny,ndump,fname)

us=zeros((Nx+1)*(Ny+1)*2,ndump);

for idump=1:ndump
    disp(['load ',num2str(idump)]);
%   dstr=sprintf('r5500/%06d.nx63.ny63.dat',idump);
    dstr=sprintf(fname,idump);
    us(:,idump)=dlmread(dstr);
end
