function invops(Bx,By,Bbpx,Bbpy)
    Bxi=kinv(Bx);
    Byi=kinv(By);
    Bbpxi=kinv(Bbpx);
    Bbpyi=kinv(Bbpy);
    return Bxi,Byi,Bbpxi,Bbpyi
end

