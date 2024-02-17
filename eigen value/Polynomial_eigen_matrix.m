
%   把一个首项系数为1的n次多项式转化为矩阵A, |lambda*E-A|=x^n+an_1*x^n-1+...+a0
%   a包含了多项式方程的系数,a0~an-1, a(i)=a_i-1
%   构造a时，要从低阶系数向高阶系数输入
function [A] = Polynomial_eigen_matrix(a)
    [n,a]=Check_column_vector(a);
    
    A=zeros(n);
    A=A+diag(repelem(1,n-1),-1);
    A(:,n)=-a(:,1);
end

