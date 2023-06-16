plot(Nsj,itsj); hold on
plot(Nsc,itsc);
plot(Nsd,itsd);

plot(Nsj,Nsj,'--');
plot(Nsj,sqrt(Nsj),'--');

xlabel('N');
ylabel('Iterations');

title('Multigrid Iterations vs. N for Spectral Discretization');
legend('Jacobi, \sigma=4/3','Chebyshev, \lambda_{min}=\lambda_{max}/4','Chebyshev, \lambda_{min}=\lambda_{max}/10','O(N)','O(N^(1/2)');
