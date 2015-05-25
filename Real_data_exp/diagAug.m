% diagonal augmentation
function [augA] = diagAug(A)
    [n, ~] = size(A);
    s = sum(A, 1);
    augA = diag(s)/(n-1) + A;
    

end