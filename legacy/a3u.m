function au=a3u(az,ay,ax,u)

% fast-tensor implementation of (az (x) ay (x) ax) u

  [mx,nx]=size(ax); [my,ny]=size(ay); [mz,nz]=size(az);
  
  if (mx*mz>=nx*nz)
    au=a3u_c2f(ax,mx,nx,ay,my,ny,az,mz,nz,u);
  else
    au=a3u_f2c(ax,mx,nx,ay,my,ny,az,mz,nz,u);
  end
end

function au=a3u_c2f(ax,mx,nx,ay,my,ny,az,mz,nz,u) % coarse-to-fine
  t1=reshape(u,nx*ny,nz); t2=zeros(nx*my,nz);
  for i=1:nz; t2(:,i)=reshape(reshape(t1(:,i),nx,ny)*ay',nx*my,1); end
  t3=reshape(t2,nx*my,nz)*az';
  t4=ax*reshape(t3,nx,my*mz);
  au=reshape(t4,mx*my*mz,1);
end

function au=a3u_f2c(ax,mx,nx,ay,my,ny,az,mz,nz,u) % fine-to-coarse
  t1=ax*reshape(u,nx,ny*nz);
  t2=reshape(t1,mx*ny,nz)*az';
  t3=zeros(mx*my,mz);
  for i=1:mz; t3(:,i)=reshape(reshape(t2(:,i),mx,ny)*ay',mx*my,1); end
  au=reshape(t3,mx*my*mz,1);
end
