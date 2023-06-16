function [p1,p2,pm,pv,psf,pp,pdiv,plot_m1,plot_m2,plot_mm,plot_m3]=setplots(Xn,Yn,Xn2,Yn2,Xm,Ym,Xpn,Ypn,Jxnl,Jynl,Jxn2l,Jyn2l,Jxml,Jyml,Jpxl,Jpyl,Xl,Yl,pdim,vort,sf,Db)

n2=round(length(reshape(Xn,[],1)));

plot_m1=@(u) mplot_ref(u,Xn,Yn,Jxnl,Jynl,Xl,Yl,pdim);
plot_mm=@(u) mplot_ref(u,Xm,Ym,Jxml,Jyml,Xl,Yl,pdim);
plot_m2=@(u) mplot_ref(u,Xpn,Ypn,Jpxl,Jpyl,Xl,Yl,pdim);
plot_m3=@(u) mplot_ref(u,Xn2,Yn2,Jxn2l,Jyn2l,Xl,Yl,pdim);

p1=@(uv)  plot_m1(uv(1:n2));
p2=@(uv)  plot_m1(uv(n2+1:end));
pm=@(uv)  plot_m1(sqrt(uv(1:n2).^2+uv(n2+1:end).^2));
pv=@(uv)  plot_m1(vort(uv));
psf=@(uv,ubsf) plot_m1(sf(uv,ubsf));

pp=@(p)    plot_m2(p);
pdiv=@(uv) plot_m2(Db(uv));
