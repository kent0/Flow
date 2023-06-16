function [uvtxy,time,dt,nu,nx,ny,ndim,nps,ngeo]=bread(fname,ifxy)

fid=fopen(fname,'r');
%[time dt nu]=fread(fid,[1 3],'double');
tmp123=fread(fid,[1 3],'double');
time=tmp123(1);
dt=tmp123(2);
nu=tmp123(3);

%[nx ny ndim nps ngeo]=fread(fid,[1 5]);
tmp12345=fread(fid,[1 5]);
nx=tmp12345(1);
ny=tmp12345(2);
ndim=tmp12345(3);
nps=tmp12345(4);
ngeo=tmp12345(5);

uv=fread(fid,[(nx+1)*(ny+1) ndim],'double');
ps=[];
xy=[];
if nps > 0; t=fread(fid,[(nx+1)*(ny+1) nps],'double'); end
if ngeo > 0; xy=fread(fid,[(nx+1)*(ny+1) ngeo],'double'); end
fclose(fid);

if ngeo == 0
    if nps == 0; uvtxy=uv; else; uvtxy=[uv,t]; end
else
    if nps == 0; uvtxy=[uv,xy]; else; uvtxy=[uv,t,xy]; end
end
