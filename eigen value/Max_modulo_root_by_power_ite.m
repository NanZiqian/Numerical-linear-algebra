
%MAXMODULOROOT_POWERITE 向量a包含了多项式方程的系数,a0-an-1, a(i)=ai-1
function [lambda1,lambda2] = Max_modulo_root_by_power_ite(a)
    a=Check_column_vector(a);
    n=size(a,1);

    A=Polynomial_eigen_matrix(a);
    [lambda1,~,flag,lambda2]= Max_modulo_eigen_by_power_ite(A);
end

