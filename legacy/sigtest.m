%figure; plot(iits1); hold on; plot(iits2); plot(iits);
%xlabel('Newton Iteration'); ylabel('GMRES Iterations');
%legend('Re=500','Re=1000','Re=1000, variable \sigma');
%
%figure; semilogy(ress1); hold on; semilogy(ress2); semilogy(ress);
%xlabel('Newton Iteration'); ylabel('Residual');
%legend('Re=500','Re=1000','Re=1000, variable \sigma');

%figure; hold on; plot(iits2); plot(iits1);
%xlabel('Newton Iteration'); ylabel('GMRES Iterations');
%legend('standard','w/ projection');

%figure; semilogy(ress2); hold on; semilogy(ress1,'--');
%xlabel('Newton Iteration'); ylabel('Residual');
%legend('standard', 'w/ projection');

%figure; hold on; plot(iits2); plot(iits1); plot(iits)
%xlabel('Newton Iteration'); ylabel('GMRES Iterations');
%legend('standard','w/ projection','w/ prec');

%figure; semilogy(ress2); hold on; semilogy(ress1,'--'); semilogy(ress);
%xlabel('Newton Iteration'); ylabel('Residual');
%legend('standard', 'w/ projection','w/ prec');

figure; plot(iits1); hold on; plot(iitsp); plot(iitsc);
xlabel('Newton Iteration'); ylabel('GMRES Iterations');
legend('w/o prec','w/ poly prec','w/ Cheby prec');

figure; semilogy(ress1); hold on; semilogy(ressp); semilogy(ressc);
xlabel('Newton Iteration'); ylabel('Residual');
legend('w/o prec', 'w/ poly prec','w/ Cheby prec');
