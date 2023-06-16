function dp=ssolve(bdti1,Bxi,Byi,Bbpxi,Bbpyi,DT,D,res)
    s=@(x) bdti1*D(a2u(Byi,Bxi,DT(x)));
    M=@(x) a2u(Bbpyi,Bbpxi,x);
    dp=pcg(s,res,1.e-6,10000,M);
end
