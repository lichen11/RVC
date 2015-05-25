function [embed] = svdembed(d, A)
% Input: test index, dimension, graph data, tau
%tau = double(tau);
%[ntotal,~] = size(A);
%ntotal = 1:ntotal;
%train = ntotal(setdiff(1:length(ntotal),test));

[U, S, V] = svd(A);
V = V(:,1:d);
Ssqrt = sqrt(S(1:d,1:d));
embed = V*Ssqrt;

%end

