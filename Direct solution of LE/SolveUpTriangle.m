function [b] = SolveUpTriangle(U,b)
% UX=b where U is a upper triangle matrix% input U and b
% output x
    N = size(U);
    n = N(1);
    for j = n:-1:2
        b(j) = b(j)/U(j,j);
        b(1:j-1)=b(1:j-1)-b(j)*U(1:j-1,j);
    end
    b(1)=b(1)/U(1,1);
