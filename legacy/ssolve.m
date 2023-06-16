function [p,flag,relres,iter,resvec]=ssolve(D,Ainv,DT,res,M,tol)
    S=@(x) D(Ainv(DT(x)));
    [p,flag,relres,iter,resvec]=pcg(S,res,tol,2000,M,[],[]);
end
