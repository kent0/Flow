function [x,rs] = cheby2(A, b, x0, iterNum, lMax, lMin)

  rs=zeros(iterNum+1,1);
  d = (lMax + lMin) / 2;
  c = (lMax - lMin) / 2;
% preCond = eye(size(A)); % Preconditioner
  x = x0;
  r = b - A * x;
  rs(1)=norm(r);

  for i = 1:iterNum % size(A, 1)
%     z = linsolve(preCond, r);
      z=r;
      if (i == 1)
          p = z;
          alpha = 1/d;
      else if (i == 2)
          beta = (1/2) * (c * alpha)^2;
          alpha = 1/(d - beta / alpha);
          p = z + beta * p;
      else
          beta = (c * alpha / 2)^2;
          alpha = 1/(d - beta / alpha);
          p = z + beta * p;
      end;

      x = x + alpha * p;
      r = b - A * x; %(= r - alpha * A * p)
      rs(i+1)=norm(r);
      if (norm(r) < 1e-15), break; end; % stop if necessary
  end;
end
