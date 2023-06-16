logk=6;
%logk=5;

%itsj=zeros(logk,1);
%Nsj=zeros(logk,1);
%stype=1;
%
%for ii=1:logk
%    N=2^(ii)+1;
%    N-1
%    mgrid1
%    itsj(ii)=it;
%    Nsj(ii)=N-1;
%end
%
%itsc=zeros(logk,1);
%Nsc=zeros(logk,1);
%stype=2;
%efac=1/4;
%
%for ii=1:logk
%    N=2^(ii)+1;
%    N-1
%    mgrid1
%    itsc(ii)=it;
%    Nsc(ii)=N-1;
%end
%
%itsd=zeros(logk,1);
%Nsd=zeros(logk,1);
%stype=2;
%efac=1/10;
%
%for ii=1:logk
%    N=2^(ii)+1;
%    N-1
%    mgrid1
%    itsd(ii)=it;
%    Nsd(ii)=N-1;
%end

%itse=zeros(logk,1);
%Nse=zeros(logk,1);
%stype=2;
%efac=1/20;
%
%for ii=1:logk
%    N=2^(ii)+1;
%    N-1
%    mgrid1
%    itse(ii)=it;
%    Nse(ii)=N-1;
%end

plot(Nsj,itsj); hold on
plot(Nsc,itsc);
plot(Nsd,itsd);
plot(Nse,itse);

plot(Nsj,Nsj,'--');
plot(Nsj,sqrt(Nsj),'--');

xlabel('N');
ylabel('Iterations');

title('Multigrid Iterations vs. N for Spectral Discretization');
%legend('Jacobi, \sigma=4/3','Chebyshev, \lambda_{min}=\lambda_{max}/4','Chebyshev, \lambda_{min}=\lambda_{max}/10');
legend('Jacobi, \sigma=4/3','Chebyshev, \lambda_{min}=\lambda_{max}/4','Chebyshev, \lambda_{min}=\lambda_{max}/10','Chebyshev, \lambda_{min}=\lambda_{max}/20','O(N)','O(N^{1/2})');


