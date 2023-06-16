function ue=ldc_reg_uv(x,y,N)

ue=[heaviside(y-1.e-4).*(x.*(1-x)).^2;y*0];
