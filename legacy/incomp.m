function [uf,p]=incomp(usf,ps,bdti1,Binv,Einv,Db,DT,RT)

dp=-ortho(Einv(Db(usf))*bdti1);
p=ps+dp;
uf=usf+RT(Binv(DT(dp))*(1./bdti1));
