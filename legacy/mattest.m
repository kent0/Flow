%ea=eig(amat);
%eb=eig(bmat);
%ep=eig(pmat);
%ec1=eig(c1mat);
%ec2=eig(c2mat);
%acmat=bmat*(amat+c1mat+c2mat);
%jmat=pmat*acmat;
%ej=eig(jmat);

%hmmat=jmat*inv(acmat);
%ehm=eig(hmmat);

%mhmat=inv(acmat)*jmat;
%emh=eig(hmmat);

%figure;
%scatter(real(ej),imag(ej))
%xlabel('Re');ylabel('Im');

%figure;
%scatter(real(emh),imag(emh));
%xlabel('Re');ylabel('Im');

%figure;
%scatter(real(ep),imag(ep));
%xlabel('Re');ylabel('Im');

%m1=jmat*inv(bmat*(amat+c2mat));
%e1=eig(m1);
%figure; scatter(real(e1),imag(e1));
%
%m2=amat+c2mat;
%e2=eig(m2);
%figure; scatter(real(e2),imag(e2));
%
%m3=bmat*(amat+c2mat);
%e3=eig(m3);
%figure; scatter(real(e3),imag(e3));

m4=amat;
e4=eig(m4);
figure; scatter(real(e4),imag(e4));
