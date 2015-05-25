function x = modifiedOMP(y, Phi, K)
% Othogonal matching pursuit of Tropp et Gilbert
% y : data
% Phi : sensing matrix
% K : sparsity

[~,N] = size(Phi);
x = zeros(1,N);
S = [];         % positions indexes of components of s
res = y;        % first residual
PhiS = [];      % Matrix of the columns used to represent y
PhiSquare = Phi^2;
for t=1:K;
    [~,j]=max(abs((PhiSquare)'*res));
    S = [S j];
    PhiS = [PhiS Phi(:,j)];
    x_est = pinv(PhiS)*y;
    res = y- PhiS*x_est;
    x(S) = x_est;
end;



end

