function mplot_ref(un,xn,yn,Jx,Jy,xm,ym,pdim);

[mx,nx]=size(Jx);
[my,ny]=size(Jy);

um=a2u(Jy,Jx,un(1:(nx*ny),1));
um=reshape(um,mx,my);
vn=reshape(un(1:(nx*ny),1),nx,ny);

if pdim == 2
    clevels=[linspace(min(un)*1.001,-1.e-4,10) -1.e-5 -1.e-6 -1.e-7 -1.e-8 1.e-12 3.e-12 1.e-11 3.e-11 1.e-10 3.e-10 1.e-9 3.e-8 1.e-8 3.e-8 1.e-7 3.e-7 1.e-6 3.e-6 1.e-5 3.e-5 1.e-4 3.e-4 1.e-3 3.e-3 1.e-2];
%   contourf(xm,ym,um,round(mx*1.25),'edgecolor','none');
%   fh = figure('Menu','none','ToolBar','none');
    colormap(jet(round(mx*1.5)))
    contourf(xm,ym,um,round(mx*1),'edgecolor','none');
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    axis equal;
    set(gca,'LooseInset',get(gca,'TightInset'));
    set(gcf,'Position',[100 100 500 500])
%   contourf(xm,ym,um,clevels);
else

    surf(xm,ym,um,'edgecolor','none');

    hold on;
    Ix=speye(size(Jx,2));
    Iy=speye(size(Jy,2));

    umn=a2u(Iy,Jx,un); umn=reshape(umn,mx,ny);
    X=reshape(xn,[],1);
    Y=reshape(yn,[],1);

    umn=Jx*reshape(un,nx,ny);
    xmn=Jx*xn;
    ymn=Jx*yn;

    unm=reshape(un,nx,ny)*Jy';
    xnm=xn*Jy';
    ynm=yn*Jy';

    plot3(xnm',ynm',unm','k');
    plot3(xmn,ymn,umn','k');
    hold off;
    axis equal
end

%colormap(jet(round(mx*1.5)))
%colormap(jet(400))
%colormap(hot(400))
%colormap(gray(400))
%colorbar
