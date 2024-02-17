
%   Tk-muk*I=QkRk Tk+1 = RkQk+muk*I 
function [Q,T] = wilkinson_shift_symmetric_QR(T)
    [~,n]=size(T);
    Q=eye(n);

    d=(T(n-1,n-1)-T(n,n))/2;
    if d ~= 0
        mu=T(n,n)-T(n,n-1)^2/( d+sign(d)*sqrt( d^2+T(n,n-1)^2 ) );
    else
        mu=T(n,n);
    end
    x=T(1,1)-mu;
    z=T(2,1);
    for k=1:n-1
        [~,~,Gk]=Givens(x,z,k,k+1,n);
        T=Gk*T*Gk;
        Q=Q*Gk;
        if k < n-1
            x=T(k+1,k);
            z=T(k+2,k);
        end
    end
end

