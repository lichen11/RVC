function [ A ] = rgSample( P )
    [n, ~] = size(P);
    a = unifrnd(0, 1, 1, n*(n-1)/2);
    U = triu(ones(n), 1);
    U(U==1) = a;
    U = U + U';
    A = (U < P) + 0;
    A(logical(eye(size(A)))) = 0;
    
   
end

