function bwrite(fname,uv,t,xy,time,dt,nu)

nx=size(uv,1)-1;
ny=size(uv,2)-1;
ndim=size(uv,3);

nps=round(length(t)/((nx+1)*(ny+1)));
ngeo=round(length(xy)/((nx+1)*(ny+1)));

fid=fopen(fname,'w');
fwrite(fid,[time dt nu],'double');
fwrite(fid,[nx ny ndim nps ngeo]);
fwrite(fid,reshape(uv,(nx+1)*(ny+1),ndim),'double');
if nps > 0; fwrite(fid,reshape(uv,(nx+1)*(ny+1),nps),'double'); end
if ngeo > 0; fwrite(fid,reshape(xy,(nx+1)*(ny+1),ngeo),'double'); end
fclose(fid);
