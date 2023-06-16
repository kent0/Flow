re=round(1./nu);
fname=['r',num2str(re),'.n',num2str(Nx),'.dat'];
fid=fopen(fname,'w'); fprintf(fid,'%20.16e\n',reshape(uf,[],1)); fclose(fid);
