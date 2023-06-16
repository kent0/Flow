function [a,b,c,x0]=ellipse(es);

m=length(es);
rmax=max(real(es));
rmin=min(real(es));

a=(rmax-rmin)*.5;
a2=a*a;
x0=(rmax+rmin)*.5;

b2=0;
for i=1:m
    t=(real(es(i))-x0)^2/a2;
    if t ~= 1
        b2=max(b2,imag(es(i))^2/(1-t));
    end
end

b=sqrt(b2);
c=sqrt(a*a-b*b);

%t=linspace(0,2*pi,1000);
%scatter(a*cos(t)+x0,b*sin(t));
