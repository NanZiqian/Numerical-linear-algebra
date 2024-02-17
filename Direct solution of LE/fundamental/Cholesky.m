
%A=LLT
function A = Cholesky(A)
    N=size(A);
    n=N(1);
    for k=1:n
        A(k,k)=sqrt(A(k,k));
        A(k+1:n,k)=A(k+1:n,k)/A(k,k);
        for j=k+1:n
            A(j:n,j)=A(j:n,j)-A(j:n,k)*A(j,k);
        end
    end
