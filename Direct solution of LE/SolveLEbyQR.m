
%用QR分解解线性方程组
function [x] = SolveLEbyQR(A,b)
       [Q,R]=QRbyHouseholder(A);
       b=Q*b;
       x=SolveUpTriangle(R,b);
end

