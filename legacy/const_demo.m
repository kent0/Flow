uk=p2k(snapsd,za,Ab);
umax=max(uk,[],2);
umin=min(uk,[],2);

delta=.05;
%uk=uk-(umax+umin)*.5;
%uk=((1-delta)*2./(umax-umin)).*uk;

%close all
%hold on
%
%plot(uk(1,:),uk(2,:));
%axis equal
%
%nk=size(uk,2);
%f=@(x,n) sum(abs(x).^n,1)-1;
%f2=@(x,y,n) (abs(x).^n+abs(y).^n)-1;
%f3=@(x,y,z,n) (abs(x).^n+abs(y).^n+abs(z).^n)-1;
%
%n=1000;
%
%fimplicit(@(x,y) f2(x,y,n));
%
%while max(f(uk(1:2,:),n*.99)) < 0; n=n*.99; end
%disp(n);
%
%fimplicit(@(x,y) f2(x,y,n));
%
%figure
%hold on
%
%n=1000;
%
%fimplicit3(@(x,y,z) f([x,y,z]',n)');
%
%while max(f(uk(18:20,:),n*.99)) < 0; n=n*.99; end
%disp(n);
%
%fimplicit3(@(x,y,z) f([x,y,z]',n)','EdgeColor','none','FaceAlpha',.5);
%%plot3(uk(1,:),uk(2,:),uk(3,:));
%plot3(uk(18,:),uk(19,:),uk(20,:));
%
%while max(f(uk(15:16,:),n*.99)) < 0; n=n*.99; end
%disp(n);
