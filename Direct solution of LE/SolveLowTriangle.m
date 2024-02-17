
%LX=b where L is a lower triangle matrix% input L and b
% output x
function [b] = SolveLowTriangle(L,b)
    [~,n]=size(L);
    for j = 1:n-1
        b(j) = b(j)/L(j,j);
        b(j+1:n)=b(j+1:n)-b(j)*L(j+1:n,j);
    end
        b(n)=b(n)/L(n,n);
