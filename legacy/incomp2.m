function uf=incomp2(usf,Binv,Einv,Db,DT,RT)

uf=usf+RT(Binv(DT(-ortho(Einv(Db(usf))))));
