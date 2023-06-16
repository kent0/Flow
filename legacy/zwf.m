      function [z,w] = zwf(N);

%     computes the N uniform nodes z and weights on [-1,1)

z=(0:N-1)'; z = -1 + 2*z./N;
w=ones(N,1); w=2*w./N;
