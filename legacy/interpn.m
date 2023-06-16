function J1=interpn(n)

xo=zwgll(n+1);
xi=zwgll(n/2+1);
J=interp_mat(xo,xi);
J1=J(2:end-1,2:end-1);
