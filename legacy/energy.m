function en=energy(x,y,H,b);
    en=.5*(x.*x*H(1,1)+y.*y*H(2,2)+x.*y*(2*H(1,2)))-b(1)*x-b(2)*y;
