function [coefficients, recoveredMat] = srclust(P, s)
    [n, ~] = size(P);
    Phi = bsxfun(@times,P, 1./sqrt(sum(P.^2, 1))); 
    %Phi = P;
    Phi(isnan(Phi)) = 0 ;
    coefficients = zeros(n, s);
    recoveredMat = zeros(n-1, n);
    for i = 1:n
        X = OMP(P((1:n) ~= i, i), Phi((1:n) ~= i, (1:n) ~= i), s);
        coefficients(i, :) = X(X ~= 0);
        selected_vertices = (X ~= 0);
        recoveredMat(:, i) = Phi(selected_vertices, (1:n) ~= i)*X';
        
    
    end
    
    
end

