n=round(size(cs,1)*.5);
cs1=cs(:,100);
nc=nb+1;
[PUi,ps,U]=deim(cb_vec(1:n2,:),nc);
c=PUi*cs1(ps');
r=cs1-cb_vec(:,1:nc)*c;
wf=cs1(:,1)*0;
wf(ps')=1.;

figure; plot_m3(wf); xlabel('x'); ylabel('y'); title('Interpolation Points');

%figure; plot_m3(sqrt(cb_vec(1:n2,1).^2+cb_vec(n2+1:end,1).^2)); xlabel('x'); ylabel('y'); title('Mode 1 Magnitude');
%figure; plot_m3(sqrt(cb_vec(1:n2,2).^2+cb_vec(n2+1:end,2).^2)); xlabel('x'); ylabel('y'); title('Mode 2 Magnitude');
%figure; plot_m3(sqrt(cb_vec(1:n2,4).^2+cb_vec(n2+1:end,4).^2)); xlabel('x'); ylabel('y'); title('Mode 4 Magnitude');
%figure; plot_m3(sqrt(cb_vec(1:n2,6).^2+cb_vec(n2+1:end,6).^2)); xlabel('x'); ylabel('y'); title('Mode 6 Magnitude');

nn=nb*2;
nn=nb+1;
ccc=proj2(cs1,cb_vec,@(u)a2u(Bbyn2,Bbxn2,u),nn);

err=zeros(nn,1);
ef=cs1;

cl2=dot(cs1,a2u(Bbyn2,Bbxn2,cs1));

for i=1:nn
    ef=ef-ccc(i)*cb_vec(:,i);
    err(i)=sqrt(dot(ef,a2u(Bbyn2,Bbxn2,ef))/cl2);
end
