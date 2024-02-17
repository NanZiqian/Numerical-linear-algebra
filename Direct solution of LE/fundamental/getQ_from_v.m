
%Q=IQ1...Qr,向量v记录了r次的列变换，v(k)代表k列与v(k)列进行列变换。
%n is the dimension of matrix, size(v) is r-全主元循环次数
function [A] = getQ_from_v(v,n)
     A=repelem(1,n);
     A=diag(A);
     temp=size(v);
    for i=1:temp(2)
        A=column_trans(i,v(i),A);
    end