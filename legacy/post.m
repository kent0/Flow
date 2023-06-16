function [diff,err]=post(uf,ue,ud,p,ubsf,t,dt,Xn,Yn,istep,iostep,nsteps,ifexact,ifplot, ...
              Bb,p1,p2,pm,pv,psf,pp,pdiv);

diff=sqrt(dot(ud,Bb(ud))/dot(uf,Bb(uf)))/dt;
err=0;
if mod(istep,iostep)==0 || istep == nsteps
    idump=round(istep/iostep);
    cstr=sprintf('%06d',idump);
    cfl=gcfl(uf,dt,Xn,Yn);
    str=['TIME=',num2str(t,'%8.3e'),'  DT=',num2str(dt,'%8.3e'),'  CFL=',num2str(cfl,'%8.3e'),'  DIFF=',num2str(diff,'%8.3e')];

    if (ifexact)
        uerr=uf-ue;
        err=sqrt(dot(uerr,Bb(uerr))/dot(ue,Bb(ue)));
        str=[str,'  ERR=',num2str(err,'%8.3e')];
    end

    disp(str);

    if ifplot || istep == nsteps
%   if false
        hold off;
%       size(uf)
%       size(ubsf)
%       psf(uf,ubsf); title(['Stream Function  ' str]);
%       hold on; qunit(uf,Xn,Yn);
%       pv(uf); title(['Vorticity  ' str]);

        pm(uf); %title(['Velocity Magnitude  ' str],'FontName','Operator Mono');

%       p1(uf-ue); title(['X-Velocity  ' str]);
%       p2(uf); title(['Y-Velocity  ' str]);

%       pp(p); title(['Pressure  ' str]);
%       pdiv(uf); title(['Divergence  ' str]);

%       xlabel('x'); ylabel('y');
%       set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12.8 7.2])

%       savefig(['plot',cstr,'.fig']);
%       saveas(gcf,['plot',cstr],'png');
%       print(['plot',cstr],'-dpng','-r100');
        saveas(gcf,['plot',cstr,'.eps'],'epsc');
%       drawnow;
    end

    [nx1,ny1]=size(Xn);
    ux=[uf;reshape(Xn,[],1);reshape(Yn,[],1)];
    ux=reshape(ux,nx1*ny1,4);
%   fid=fopen([cstr,'.nx',num2str(nx1-1),'.ny',num2str(ny1-1),'.dat'],'w'); fprintf(fid,'%20.16e %20.16e %20.16e %20.16e\n',ux); fclose(fid);
    fid=fopen([cstr,'.nx',num2str(nx1-1),'.ny',num2str(ny1-1),'.dat'],'w'); fprintf(fid,'%20.16e\n',uf); fclose(fid);
end
