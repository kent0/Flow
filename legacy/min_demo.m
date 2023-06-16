close all;
cmap = colormap(lines(7));
nclines=80;

base='/Users/kaneko/Developer/Cases/ldc_t0/run/';
uk=dlmread([base 'ops/uk']);
lb=round(dlmread([base 'ops/nb']));
H1=dlmread([base 'ops/h1']);
b1=dlmread([base 'ops/b1']);

nb=length(b1);
H1=reshape(H1,nb,nb);
[VH,DH]=eig(H1);
ns=1000;
uk=reshape(uk,lb+1,ns);

umin=min(uk,[],2);
umax=max(uk,[],2);

u1=uk(2,:);
u2=uk(3,:);

us=H1\b1;
nbox=1000;

u1min=min(u1);
u1max=max(u1);
u2min=min(u2);
u2max=max(u2);

[xx,yy]=meshgrid(linspace(-10,10,1000),linspace(-10,10,1000));

x=[linspace(u1min,u1max,nbox)';repmat(u1max,nbox,1);linspace(u1max,u1min,nbox)';repmat(u1min,nbox,1)];
y=[repmat(u2min,nbox,1);linspace(u2min,u2max,nbox)';repmat(u2max,nbox,1);linspace(u2max,u2min,nbox)'];

plot(x,y,'color',cmap(1,:)); hold on
plot(0,0,'.','color',cmap(1,:));
%plot(us(1),us(2),'x','color',cmap(1,:));

b2=b1-H1*[0;0;us(3:end)]; b2=b2(1:2);
H2=H1(1:2,1:2);
f=@(x,y) energy(x,y,H2,b2);
zz=f(xx,yy); minz=min(min(zz));
contour(xx,yy,f(xx,yy),linspace(0,10,nclines).^2+minz,'color',cmap(1,:),'linewidth',2);

%fcontour(f,[-2 8 -5 5]);
axis equal

[minbox,ind]=min(f(x,y));
plot(x(ind),y(ind),'*','color',cmap(1,:));
plot(u1,u2,'k.');

xlabel('$u_1$','Interpreter','latex');
ylabel('$u_2$','Interpreter','latex');

[W,lam]=eig(H2);
lam=flip(flip(lam,1),2);
W=flip(W,2);
W=W*(-1);
%W lam WT u = b
%lam uh = bh;
bh=W'*b2;
fh=@(x,y) energy(x,y,lam,bh);

figure; hold on
zz=fh(xx,yy); minz=min(min(zz));
contour(xx,yy,fh(xx,yy),linspace(0,10,nclines).^2+minz,'color',cmap(2,:),'linewidth',2);

uh=W'*[u1;u2];

uhmax=max(uh,[],2);
uhmin=min(uh,[],2);

xh=[linspace(uhmin(1),uhmax(1),nbox)';repmat(uhmax(1),nbox,1);linspace(uhmax(1),uhmin(1),nbox)';repmat(uhmin(1),nbox,1)];
yh=[repmat(uhmin(2),nbox,1);linspace(uhmin(2),uhmax(2),nbox)';repmat(uhmax(2),nbox,1);linspace(uhmax(2),uhmin(2),nbox)'];

[minbox,indh]=min(fh(xh,yh));
t1=W'*[x(ind);y(ind)];
plot(t1(1),t1(2),'*','color',cmap(1,:));
plot(xh(indh),yh(indh),'*','color',cmap(2,:));

plot(xh,yh,'color',cmap(2,:));

tmp=W'*[x';y']; plot(tmp(1,:),tmp(2,:),'color',cmap(1,:));
axis equal

plot(uh(1,:),uh(2,:),'k.');

xlabel('$\hat{u}_1$','Interpreter','latex');
ylabel('$\hat{u}_2$','Interpreter','latex');

%W lam WT u = b
%lam uh = bh;

bc=(lam^(-.5))*bh;
fc=@(x,y) energy(x,y,eye(2),bc);

figure; hold on
zz=fc(xx,yy); minz=min(min(zz));
contour(xx,yy,fc(xx,yy),linspace(0,10,nclines).^2+minz,'color',cmap(3,:),'linewidth',2);

uc=lam^(.5)*uh;

ucmax=max(uc,[],2);
ucmin=min(uc,[],2);

xc=[linspace(ucmin(1),ucmax(1),nbox)';repmat(ucmax(1),nbox,1);linspace(ucmax(1),ucmin(1),nbox)';repmat(ucmin(1),nbox,1)];
yc=[repmat(ucmin(2),nbox,1);linspace(ucmin(2),ucmax(2),nbox)';repmat(ucmax(2),nbox,1);linspace(ucmax(2),ucmin(2),nbox)'];

tmp=(lam^.5)*W'*[x';y']; plot(tmp(1,:),tmp(2,:),'color',cmap(1,:));
tmp=(lam^.5)*[xh';yh']; plot(tmp(1,:),tmp(2,:),'color',cmap(2,:));
plot(xc,yc,'--','color',cmap(3,:));

t1=(lam^.5)*W'*[x(ind);y(ind)];
plot(t1(1),t1(2),'*','color',cmap(1,:));
t1=(lam^.5)*[xh(indh);yh(indh)];
plot(t1(1),t1(2),'*','color',cmap(2,:));
[minbox,indc]=min(fc(xc,yc));
plot(xc(indc),yc(indc),'*','color',cmap(3,:));

plot(uc(1,:),uc(2,:),'k.');

xlabel('$\check{u}_1$','Interpreter','latex'); ylabel('$\check{u}_2$','Interpreter','latex');

axis equal
