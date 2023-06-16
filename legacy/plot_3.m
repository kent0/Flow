        figure
        m=length(u2);
        time=(1:m)*100/m;
        semilogy(time,u2); hold on
        semilogy(time,uinf);
        semilogy(time,up2);
        semilogy(time,upinf);
        xlabel('Time'); ylabel('||u||'); hold on
        legend("||u||_2","||u||_\infty","||u'||_2","||u'||_\infty");
