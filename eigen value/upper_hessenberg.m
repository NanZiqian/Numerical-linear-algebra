
%   Q^H A Q = H
function [Q,H] = upper_hessenberg(A)
    [~,n]=size(A);
    Q=eye(n);

    for k = 1:n-2
        [v,beta]=House(A(k+1:n,k));
        temp=(eye(n-k)-beta*v*v');
        Q=Q*blkdiag(eye(k),temp);
        A(k+1:n,k:n)=temp*A(k+1:n,k:n);
        A(1:n,k+1:n)=A(1:n,k+1:n)*temp;
    end

    H=A;
end

