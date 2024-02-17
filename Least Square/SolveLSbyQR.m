
%SOLVELSBYQR 此处显示有关此函数的摘要
function [x] = SolveLSbyQR(A,b)
    N=size(A);
    m=N(1);
    n=N(2);
    %% QR分解A
    [Q,R]=QRbyHouseholder(A);
    Q1=Q(:,1:n);
    %% 计算c1
    c1=Q1'*b;
    %% 求解Rx=c1
    R1=R(1:n,:);
    x=SolveUpTriangle(R1,c1);
end

