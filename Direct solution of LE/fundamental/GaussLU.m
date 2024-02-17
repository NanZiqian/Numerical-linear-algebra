function A = GaussLU(A,N)
%对A通过Gauss消元LU分解N= size(A);
    n = N(1);
    for k = 1:n-1
        A(k+1:n,k) = A(k+1:n,k)/A(k,k);
        A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-A(k+1:n,k)*A(k,k+1:n);
    end
