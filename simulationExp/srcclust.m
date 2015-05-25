function [coefficients, residual] = srcclust(P, s)
    [n, ~] = size(P);
    Phi = bsxfun(@times,P, 1./sqrt(sum(P.^2, 1))); 
    %Phi = P;
    Phi(isnan(Phi)) = 0 ;
    coefficients = zeros(n, s);
    residual = zeros(n,1);
    %selected_verticesMat = cell(n, 1);
    %recoveredMat = zeros(n-1, n);
    %h = waitbar(0,'Please wait...');
    for i = 1:n
    %    waitbar(i/n)
        X = modifiedOMP(P((1:n) ~= i, i), Phi((1:n) ~= i, (1:n) ~= i), s);
        coefficients(i, :) = X(X ~= 0);
        residual(i) = norm(Phi((1:n) ~= i, (1:n) ~= i) * X' - P((1:n) ~= i, i), 2);
        %selected_verticesMat{i} = (X ~= 0);
        %recoveredMat(:, i) = Phi(selected_vertices, (1:n) ~= i)*X';
        
    
    end
    %close(h);
    
end

