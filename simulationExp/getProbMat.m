function [P, tau] = getProbMat(B, rho, n)
k = length(rho);
tau =  randsample(k, n, 1, rho);
[vec, val] = eigs(B, k);
centers = vec * sqrt(val);
latent = centers(tau,:);
P = latent * (latent)';
end

