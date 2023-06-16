function [f1,f2]=evalf(uf,C,Hb,Binv,Einv,Db,DT,R,RT)

ub=uf-RT(R(uf));

f1=R(incomp2(RT(Binv(R(Hb(uf,0)+C(uf,uf)))),Binv,Einv,Db,DT,RT));
f2=Binv(DT(Einv(Db(ub))));
