function [betas,alphas1,alphas2]=bdfext_var(ts,k)

alphas1=interp_mat(ts(1),ts(2:end))';
alphas2=interp_mat(ts(1),ts(2:end-1))';
betas=deriv_mat(ts)'; betas=betas(:,1);

alphas1=[alphas1; zeros(k-length(alphas1),1)];
alphas2=[alphas2; zeros(k-length(alphas2),1)];
betas=[betas; zeros(k+1-length(betas),1)];
