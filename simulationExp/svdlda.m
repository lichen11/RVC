%function [pred] = svdlda(test, d, A, tau)
function [pred] = svdlda(test, embed, tau)
% Input: test index, dimension, graph data, tau
tau = double(tau);
[ntotal,~] = size(tau);
ntotal = 1:ntotal;
train = ntotal(setdiff(1:length(ntotal),test));
%[U, S, ~] = svd(A);
%U = U(:,1:d);
%Ssqrt = sqrt(S(1:d,1:d));
%embed = U*Ssqrt;
pred = classify(embed(test,:),embed(train,:), tau(train),'linear');


