function [betas,alphas]=bdfext(k)

if (k == 0)
    betas=[0,0,0,0];
    alphas=[0,0,0];
elseif (k == 1)
    betas=[1,-1,0,0];
    alphas=[1,0,0];
elseif (k == 2)
    betas=[3/2,-2,1/2,0];
    alphas=[2,-1,0];
elseif (k == 3)
    betas=[11/6,-3,3/2,-1/3];
    alphas=[3,-3,1];
else
    disp('unsupported k in bdfext');
end
