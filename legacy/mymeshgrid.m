function [Xn,Yn,Zn]=mymeshgrid(xm,ym,zm,ndim)

switch ndim
    case 1
        Xn=xm; Yn=Xn*0; Zn=Xn*0;
    case 2
        [Xn,Yn]=meshgrid(xn,yn); Zn=Xn*0;
    case 3
        [Xn,Yn,Zn]=meshgrid(xn,yn,zn);
    otherwise
        Xn=[];Yn=[];Zn=[];
end
