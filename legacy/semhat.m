      function[Ah,Bh,Ch,Dh,z,w] =  semhat(N)
%
%                                                 ^
%     Compute the single element 1D SEM Stiffness Mass, and Convection
%     matrices, as well as the points and weights for a polynomial
%     of degree N
%

      [z,w] = zwgll(N);

      Bh    = diag(w);   Bh=sparse(Bh);
      Dh    = deriv_mat(z);

      Ah    = Dh'*Bh*Dh; Ah=0.5*(Ah+Ah');
      Ch    = Bh*Dh;

