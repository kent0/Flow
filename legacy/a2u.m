function au=a2u(ay,ax,u)

% fast-tensor implementation of (ay (x) ax) u

  [mx,nx]=size(ax); [my,ny]=size(ay);
  n=nx*ny;
  m=mx*my;
  dim=round(size(u,1)/n);
  col=size(u,2);
  au=zeros(m*dim,col);
  for ic=1:col
  for id=1:dim
    au((1+(id-1)*m):(id*m),ic)=reshape(ax*reshape(u((1+(id-1)*n):(id*n),ic),nx,ny)*ay',m,1);
  end
  end
% dim=round(length(u)/(nx*ny));
% if (dim == 1)
%   au=reshape(ax*reshape(u,nx,ny*dim)*ay',mx*my,1);
% else
%   v=reshape(u,nx*ny,2);
%   au=[reshape(ax*reshape(v(:,1),nx,ny)*ay',mx*my,1); ...
%       reshape(ax*reshape(v(:,2),nx,ny)*ay',mx*my,1)];
% end
end
