function vort=v2d(uv,Dbx,Dby)
    n2=round(length(uv)/2);
    Ibx=speye(size(Dbx,1)); Iby=speye(size(Dby,1));
    vort=a2u(Iby,Dbx,uv(n2+1:end))-a2u(Dby,Ibx,uv(1:n2));
    size(vort)
end
