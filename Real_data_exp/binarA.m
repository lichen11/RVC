function [A] = binarA(A)
    A(logical(eye(size(A)))) = 0;
    A=A~=0;
    A = double(A);
end


