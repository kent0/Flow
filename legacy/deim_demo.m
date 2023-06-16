n=round(size(cs,1)*.5);
cs1=cs(:,100);
nc=nb+1;
[PUi,ps,U]=deim(cb_vec(1:n2,:),nc);
%[PUi,ps,U]=deim(cb(1:n2,:),nb+1);
c=PUi*cs1(ps');
r=cs1-cb_vec(:,1:nc)*c;

figure; plot_m3(sqrt(r(1:n2).^2+r(n2+1:end).^2)); xlabel('x'); ylabel('y'); title('Interpolation Error Magnitude of Convective Term');
figure; plot_m3(sqrt(cs1(1:n2).^2+cs1(n2+1:end).^2)); xlabel('x'); ylabel('y'); title('Magnitude of Convective Term');

figure; ifexactc=false; rom; wf=vf;
        ifexactc=true;  rom;

xlabel('x'); ylabel('y'); title('Standard ROM Velocity Magnitude at Time T=100');

figure;
pm(wf-vf);
xlabel('x'); ylabel('y'); title('DEIM ROM Velocity Difference Magnitude at Time T=100');

