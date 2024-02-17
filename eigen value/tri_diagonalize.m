
%   只对对称矩阵有用,T=Q^H*A*Q，没有置零，因此要取三个对角的元素
function [Q,T] = tri_diagonalize(A)
    [~,n]=size(A);
    Q=eye(n);
    for k=1:n-2
        [v,beta,Hk]=House(A(k+1:n,k));
        Q=Q*blkdiag(eye(k),Hk);
        u=beta*A(k+1:n,k+1:n)*v;
        w=u-beta/2*u'*v*v;
        A(k+1,k)=norm(A(k+1:n,k),2);
        A(k,k+1)=A(k+1,k);
        A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-v*w'-w*v';
    end
    T=triu(tril(A,1),-1);
end

