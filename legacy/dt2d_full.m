function dt=dt2d_full(Bbx,Bby,Dbx,Dby,Jpx,Jpy)
    dt=[kron(Bby*Jpy,Dbx'*Bbx*Jpx);kron(Dby'*Bby*Jpy,Bbx*Jpx)];
end
