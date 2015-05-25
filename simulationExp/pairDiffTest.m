function [pvalue] = pairDiffTest(result1, result2)
%   paired difference test using Wilcoxon rank sum test
%   Input: matrices the rows are varying according to some parameters.
%          the columns are the L simulations
%   suppose we have M parameters, then we will have M p-values
%   result1 and result2 are M by L
%   I only want one-sided so I mult the pvalues by 2 in the end
[M, ~] = size(result1);
pvalue = zeros(M, 1);
for k = 1:M
  % pvalue(k) = ranksum(result1(k, : ), result2(k, :))/2;
   pvalue(k) = signrank(result1(k, : ), result2(k, :))*2;
end

end

