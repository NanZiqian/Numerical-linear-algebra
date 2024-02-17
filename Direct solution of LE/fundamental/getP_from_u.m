
% P=Pr...P1I,向量u记录了r次的行变换，u(k)代表k行与u(k)行进行行变换。
%n is the dimension of matrix, size(u) is r-全主元循环次数
function [A] = getP_from_u(u,n)
     A=repelem(1,n);
     A=diag(A);
     temp=size(u);
    for i=1:temp(2)
        A=row_trans(i,u(i),A);
    end


