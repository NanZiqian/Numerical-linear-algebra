
%need matrix A and b beforehand
function [x] = SolveLEbyGauss(A,b)
    %如果输入行向量，则自动转为列向量    
    N=size(b);
    if N(2) ~= 1
        b=b';
    end
    n = size(A,1);
    B = GaussLU(A,n);
    L = eye(n)+tril(B,-1);
    U = triu(B);
    y=SolveLowTriangle(L,b);
    x=SolveUpTriangle(U,y);
end
