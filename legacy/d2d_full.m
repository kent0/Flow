function d=d2d_full(Bbx,Bby,Dbx,Dby,Jpx,Jpy)
    d=[kron(Jpy'*Bby,Jpx'*Bbx*Dbx) kron(Jpy'*Bby*Dby,Jpx'*Bbx)];
end
