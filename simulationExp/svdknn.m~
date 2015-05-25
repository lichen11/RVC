%function [pred, embed ] = svdknn(test, d, A, tau, k)
function [pred] = svdknn(test, embed, tau, k)
%UNTITLED Summary of this function goes here
% Detailed explanation goes here
% Input: A, d, k neighbor, test, tau
% Output: embed, predicted value, 
tau = double(tau);
%[U, S, ~] = svd(A);
%U = U(:,1:d);
%Ssqrt = sqrt(S(1:d,1:d));
%embed = U*Ssqrt;

[ntotal,~] = size(tau);
ntotal = 1:ntotal;
train = ntotal(setdiff(1:length(ntotal),test));
pred = knnclassify(embed(test,:), embed(train,:), tau(train), k);
%end


