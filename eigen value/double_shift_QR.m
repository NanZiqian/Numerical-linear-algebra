
%   Detailed explanation goes here
function [Q,H] = double_shift_QR(H)
    [~,n]=size(H);
    
    m=n-1;
    s=H(m,m)+H(n,n);
    t=H(m,m)*H(n,n)-H(m,n)*H(n,m);
    x=H(1,1)*H(1,1)+H(1,2)*H(2,1)-s*H(1,1)+t;
    y=H(2,1)*(H(1,1)+H(2,2)-s);
    z=H(2,1)*H(3,2);
    Q=eye(n);
    for k=0:n-3
        [~,~,Pk]=House([x;y;z]);
        Q=Q*blkdiag(eye(k),Pk,eye(n-k-3));
        q=max([1,k]);
        H(1+k:k+3,q:n)=Pk*H(k+1:k+3,q:n);
        r=min([k+4,n]);
        H(1:r,k+1:k+3)=H(1:r,k+1:k+3)*Pk;
        x=H(k+2,k+1);
        y=H(k+3,k+1);
        if k < n-3
            z = H(k+4,k+1);
        end
    end
    [~,~,Pn_2]=House([x;y]);
    Q=Q*blkdiag(eye(n-2),Pn_2);
    H(n-1:n,n-2:n)=Pn_2*H(n-1:n,n-2:n);
    H(1:n,n-1:n)=H(1:n,n-1:n)*Pn_2;
end

