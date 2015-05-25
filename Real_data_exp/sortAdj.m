function [sortedA, sortedLab ] = sortAdj(originalA, originalLab)
%Input:
%Output:

sortedA = [originalLab';  originalA];
[d1, ~] = size(sortedA);
[~, d2] = sort(sortedA(1,:));
sortedA = sortedA(:,d2);
sortedA = sortedA(2:d1,:);
sortedA = [originalLab sortedA];
sortedA = sortrows(sortedA, 1);
sortedA = sortedA(:, 2:d1);
sortedLab = sort(originalLab);

end

