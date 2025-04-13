let ulag,clag,plag,tlag
    global function advance(uf,ub,p,t,dt,istep,Hb,Hinv,Bb,Binv,DbT,R,RT,C)
        c=C(uf,uf);

        if (istep == 1)
            n=length(uf); np=length(p); nc=length(c);
            ulag=zeros(n,3);
            clag=zeros(nc,3);
            plag=zeros(np,3);
            tlag=zeros(1,3);
        end

    #   println("ulag $(size(ulag)), $(size(uf))")
        lag!(ulag,uf);
        lag!(clag,c);
        lag!(plag,p);
        lag!(tlag,[t]);

        ts=[t+dt;tlag[1:min(3,istep)]]
        bdti,alphas1,alphas2=bdfext_var(ts,3);
        bdti1=bdti[1];

        b=-Bb(ulag*bdti[2:end]);
        b=b-clag*alphas1;
        b=b-Hb(ub,bdti1);
        ps=plag*alphas2;
        b=R(b+DbT(ps));
        
        usf=RT(Hinv(b,bdti1,istep<4))+ub;

        return usf,ps,bdti1
    end
end

function bdfext_var(ts,k)
    alphas1=interp_mat(ts[1],ts[2:end])';
    alphas2=interp_mat(ts[1],ts[2:end-1])';
    betas=deriv_mat(ts)'; betas=betas[:,1];

    alphas1=[alphas1; zeros(k-length(alphas1),1)];
    alphas2=[alphas2; zeros(k-length(alphas2),1)];
    betas=[betas; zeros(k+1-length(betas),1)];

    return betas,alphas1,alphas2
end
