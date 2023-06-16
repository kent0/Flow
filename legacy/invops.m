function [Bxi,Byi,Bbpxi,Bbpyi]=invops(Bx,By,Bbpx,Bbpy);

Bxi=inv(Bx);
Byi=inv(By);
Bbpxi=inv(Bbpx);
Bbpyi=inv(Bbpy);
