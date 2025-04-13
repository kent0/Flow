function post(uf,ue,ud,p,ubsf,t,dt,Xn,Yn,istep,iostep,nsteps,ifexact,ipstep,
              Bb,p1,p2,pm,pv,psf,pp,pdiv,pfx="");

    diff=sqrt(dot(ud,Bb(ud))/dot(uf,Bb(uf)))/dt;
    err=0;

    if mod(istep,ipstep)==0 || istep == nsteps
        iplot=round(istep/ipstep);
        dstr=@sprintf("%06d",iplot);
        pm(uf,pfx*"o"*dstr); #title(['Velocity Magnitude  ' str],'FontName','Operator Mono');
        pp(p,pfx*"p"*dstr); #title(['Velocity Magnitude  ' str],'FontName','Operator Mono');
#           println("plot...")
#           sleep(5)
    end
    
    if mod(istep,iostep)==0 || istep == nsteps
        idump=round(istep/iostep);
        cstr=@sprintf("%06d",idump);
        cfl=gcfl(uf,dt,Xn,Yn);
        str="TIME="*@sprintf("%8.3e",t)*"  DT="*@sprintf("%8.3e",dt)*"  CFL="*@sprintf("%8.3e",cfl)*"  DIFF="*@sprintf("%8.3e",diff);

#       if (ifexact)
#           uerr=uf-ue;
#           err=sqrt(dot(uerr,Bb(uerr))/dot(ue,Bb(ue)));
#           str=str*"  ERR="*num2str(err,'%8.3e')];
#       end


        nx1,ny1=size(Xn);
#       ux=[uf;reshape(Xn,[],1);reshape(Yn,[],1)];
#       ux=reshape(ux,nx1*ny1,4);
        fname=pfx*cstr*".nx$(nx1).ny$(ny1).dat"
        write(fname,uf);
        if istep == nsteps
            fname=pfx*"geo.nx$(nx1).ny$(ny1).dat"
            write(fname,[reshape(Xn,:); reshape(Yn,:)]);
        end
#       fid=fopen([cstr*".nx$(nx1-1).ny$(ny1-1).dat"],"w"); fprintf(fid,"%20.16e\n",uf); fclose(fid);
#       write(
        println(str)
    end

    return diff,err
end

function mplot_ref(un,xn,yn,Jx,Jy,xm,ym,pdim,fname="figure",pfx="")
    mx,nx=size(Jx);
    my,ny=size(Jy);

    um=a2u(Jy,Jx,un[1:(nx*ny),1]);
    um=reshape(um,mx,my);
    vn=reshape(un[1:(nx*ny),1],nx,ny);
#   println("inside mplot_ref")

#   contour(xn,yn,un)
#   println("size un $(size(un))")
#   println("size xn $(size(xn))")
#   println("size yn $(size(yn))")
#   heatmap(reshape(un,nx,ny))
#   heatmap(xn,yn,reshape(un,nx,ny))
#   heatmap(reshape(xn,:),reshape(yn,:),reshape(un,nx,ny))
#   heatmap(xn[:,1],yn[1,:],reshape(un,nx,ny)')
#   h=heatmap(xn[:,1],yn[1,:],reshape(un,nx,ny)')

#   h=heatmap(xm[:,1],ym[1,:],reshape(um,mx,my)',aspect_ratio=1.0,c=:jet1,legend=false,axis=([], false))


#   savefig(pfx*fname*".png")

    umax = maximum(um)
    umin = minimum(um)
    
    if umax == umin
        um = um .- umin
    else
        um = (um .- umin) .* ( 1 / (umax - umin))
    end

    cmap = ColorSchemes.turbo
    plotu = reshape(um,mx,my)'
    
    ucolor = map(x->get(cmap,x), plotu[end:-1:1,:])
    save(pfx*fname*".png",ucolor)

    if pdim == 2
#       clevels=[linspace(min(un)*1.001,-1.e-4,10) -1.e-5 -1.e-6 -1.e-7 -1.e-8 1.e-12 3.e-12 1.e-11 3.e-11 1.e-10 3.e-10 1.e-9 3.e-8 1.e-8 3.e-8 1.e-7 3.e-7 1.e-6 3.e-6 1.e-5 3.e-5 1.e-4 3.e-4 1.e-3 3.e-3 1.e-2];
#   %   contourf(xm,ym,um,round(mx*1.25),'edgecolor','none');
#   %   fh = figure('Menu','none','ToolBar','none');
#       colormap(jet(round(mx*1.5)))
#       contourf(xm,ym,um,round(mx*1),'edgecolor','none');
#       set(gca,'XTick',[])
#       set(gca,'YTick',[])
#       axis equal;
#       set(gca,'LooseInset',get(gca,'TightInset'));
#       set(gcf,'Position',[100 100 500 500])
#   %   contourf(xm,ym,um,clevels);
    else

#       surf(xm,ym,um,'edgecolor','none');

#       hold on;
#       Ix=speye(size(Jx,2));
#       Iy=speye(size(Jy,2));

#       umn=a2u(Iy,Jx,un); umn=reshape(umn,mx,ny);
#       X=reshape(xn,[],1);
#       Y=reshape(yn,[],1);

#       umn=Jx*reshape(un,nx,ny);
#       xmn=Jx*xn;
#       ymn=Jx*yn;

#       unm=reshape(un,nx,ny)*Jy';
#       xnm=xn*Jy';
#       ynm=yn*Jy';

#       plot3(xnm',ynm',unm','k');
#       plot3(xmn,ymn,umn','k');
#       hold off;
#       axis equal
    end
end
