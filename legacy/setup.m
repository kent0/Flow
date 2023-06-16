function [Ax,Ay,Abx,Aby,Bx,By,Bbx,Bby,Bbxm,Bbym,Bbpxi,Bbpyi, ...
          Dbx,Dby,Dbxm,Dbym,Dbpx,Dbpy,Ibx,Iby,Ibxm,Ibym,Ibpx,Ibpy, ...
          Jxm,Jym,Jxl,Jyl,Jpx,Jpy,Jpxl,Jpyl,Px,Py,Rx,Ry, ...
          Xn,Yn,Xm,Ym,Xl,Yl,Xpn,Ypn]=setup(Nx,Ny,dom,bc,nf);

    k=2;

    ax=dom(1);bx=dom(2);ay=dom(3);by=dom(4);

    lx=bx-ax;
    ly=by-ay;

    Npx=Nx-k;
    Npy=Ny-k;

    Mx=ceil(Nx*1.5);
    My=ceil(Ny*1.5);

    Lx=255;
    Ly=255;

    [zxn,wxn]=zwgll(Nx);
    [zyn,wyn]=zwgll(Ny);

    [zxm,wxm]=zwgll(Mx);
    [zym,wym]=zwgll(My);

    [zxl,wxl]=zwgll(Lx);
    [zyl,wyl]=zwgll(Ly);

    xf=ax+lx*(zwgll(Nx-nf)+1)*.5;
    yf=ay+ly*(zwgll(Ny-nf)+1)*.5;

    [zxpn,wxpn]=zwgll(Npx);
    [zypn,wypn]=zwgll(Npy);

%   [zxpn,wxpn]=zwgl(Npx+1);
%   [zypn,wypn]=zwgl(Npy+1);

    wxn=lx*wxn*.5;
    wyn=ly*wyn*.5;
    wxm=lx*wxm*.5;
    wym=ly*wym*.5;
    wxl=lx*wxl*.5;
    wyl=ly*wyl*.5;

    wxpn=lx*wxpn*.5;
    wypn=ly*wypn*.5;

    xn=ax+lx*(zxn+1)*.5;
    yn=ay+ly*(zyn+1)*.5;
    xm=ax+lx*(zxm+1)*.5;
    ym=ay+ly*(zym+1)*.5;
    xl=ax+lx*(zxl+1)*.5;
    yl=ay+ly*(zyl+1)*.5;
    xpn=ax+lx*(zxpn+1)*.5;
    ypn=ay+ly*(zypn+1)*.5;

    [Xn,Yn]=ndgrid(xn,yn);
    [Xm,Ym]=ndgrid(xm,ym);
    [Xl,Yl]=ndgrid(xl,yl);

    [Xpn,Ypn]=ndgrid(xpn,ypn);

    if bc(1) == 'd'; i1=2; elseif bc(1) == 'n'; i1=1; end
    if bc(2) == 'd'; i2=1; elseif bc(2) == 'n'; i2=0; end
    if bc(3) == 'd'; i3=2; elseif bc(3) == 'n'; i3=1; end
    if bc(4) == 'd'; i4=1; elseif bc(4) == 'n'; i4=0; end

    Rx=speye(Nx+1); Rx=Rx(i1:end-i2,:);
    Ry=speye(Ny+1); Ry=Ry(i3:end-i4,:); % Dirichlet-Dirichlet

    Bbx=sparse(diag(wxn));
    Bby=sparse(diag(wyn));

    Bx=Rx*Bbx*Rx';
    By=Ry*Bby*Ry';

    Bbxm=sparse(diag(wxm));
    Bbym=sparse(diag(wym));

    Bbpxi=sparse(diag(1./wxpn));
    Bbpyi=sparse(diag(1./wypn));

    Dbx=deriv_mat(xn);
    Dby=deriv_mat(yn);

    Dbxm=deriv_mat(xm);
    Dbym=deriv_mat(ym);

    Dx=Dbx*Rx';
    Dy=Dby*Ry';

    Dbpx=deriv_mat(xpn);
    Dbpy=deriv_mat(ypn);

    Abx=Dbx'*Bbx*Dbx; Abx=.5*(Abx+Abx');
    Aby=Dby'*Bby*Dby; Aby=.5*(Aby+Aby');

    Ax=Rx*Abx*Rx'; Ax=.5*(Ax+Ax');
    Ay=Ry*Aby*Ry'; Ay=.5*(Ay+Ay');

    Ibx=speye(Nx+1);
    Iby=speye(Ny+1);

    Ibxm=speye(Mx+1);
    Ibym=speye(My+1);

    Ibpx=speye(Npx+1);
    Ibpy=speye(Npy+1);

    Jxm=interp_mat(xm,xn);
    Jym=interp_mat(ym,yn);

    Jxl=interp_mat(xl,xn);
    Jyl=interp_mat(yl,yn);

    Jpx=interp_mat(xn,xpn);
    Jpy=interp_mat(yn,ypn);

    Jpxl=interp_mat(xl,xpn);
    Jpyl=interp_mat(yl,ypn);

    Jxf=interp_mat(xf,xn);
    Jyf=interp_mat(yf,yn);

    Px=Jxf'*Jxf;
    Py=Jyf'*Jyf;
end
