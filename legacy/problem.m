function [ufun,f,dom,ifexact]=problem(cname,nu,Ny)

if strcmp(cname,'st2')
    ufun=@(x,y,t) st2_uv(x,y,t,w,nu);
    dom=[0,1,0,1];
    f=@(x,y) [x;y]*0;
    ifexact=true;
elseif strcmp(cname,'start')
    ufun=@(x,y,t) start_uv(x,y,t,10);
    dom=[-1,1,-1,1];
    f=@(x,y) [x*0+1.5*nu;y*0];
    ifexact=true;
    ifexact=false;
elseif strcmp(cname,'kov')
    ufun=@(x,y,t) kov_uv(x,y,nu);
    dom=[-.5,2.0,-.5,1.5];
    f=@(x,y) [x;y]*0;
    ifexact=true;
elseif strcmp(cname,'walsh')
    ufun=@(x,y,t) walsh_uv(x,y,t,nu);
    dom=[0,2*pi,0,2*pi];
%   dom=[-1,1,-1,1];
    f=@(x,y) [x;y]*0;
    ifexact=true;
elseif strcmp(cname,'ldc')
    ufun=@(x,y,t) [y.^(32);y*0];
    ufun=@(x,y,t) [heaviside(y-1+1.e-3);y*0];
    dom=[0,1,0,1];
    f=@(x,y) [x;y]*0;
    ifexact=false;
elseif strcmp(cname,'ldc_reg')
    ufun=@(x,y,t) [16*heaviside(y-1+1.e-3).*(x.*(1-x)).^2;y*0];%ldc_reg_uv(x,y,t);
%   dom=[-1,1,-1,1];
    dom=[0,1,0,1];
    f=@(x,y) [x;y]*0;
    ifexact=false;
elseif strcmp(cname,'swirl')
    ufun=@(x,y,t) [2*y.*(1-x.^2);-2*x.*(1-y.^2)];
    dom=[-1,1,-1,1];
    f=@(x,y) [4*y;-4*x]*nu;
    ifexact=false;
elseif strcmp(cname,'decay')
    lam=pi*.5;
    ufun=@(x,y,t) [cos(lam*y).*exp(-lam^2*nu*t);x*0];
    dom=[-1,1,-1,1];
    f=@(x,y) [x;x]*0;
    ifexact=true;
elseif strcmp(cname,'decay2')
    ufun=@(x,y,t) 1.5*[1-y.*y;x.*0]-start_uv(x,y,t,30);
    dom=[-1,1,-1,1];
    f=@(x,y) [x;x]*0;
    ifexact=true;
end
