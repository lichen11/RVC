function [A, tau, P] = sbm(rho, B, n )

tau = randsample(1:length(rho), n, true, rho);
P = B(tau, tau);
U = rand(n);
U = triu(U,1);
U = U + U';
A = (U<P)+0;
diagIdx = 1:n+1:n^2;
A(diagIdx)=0;
tau = double(tau);
tau = tau';
end



