
function [x,L,U] = SolveLEbyGaussinall_max(A,b)
    n = size(A,1);
    % A1=A;
%% PAQ=LU
    for k=1 : n-1
        p=k;
        q=k;
        for i=k+1 : n
            for j=k+1 : n
                if abs(A(i,j))>abs(A(p,q))
                    p=i;
                    q=j;
                end
            end
        end
        %交换k行和p行，k列和q列
        A=row_trans(k,p,A);
        A=column_trans(k,q,A);
        u(k)=p;
        v(k)=q;
        if A(k,k) ~= 0
            A(k+1:n,k)=A(k+1:n,k)/A(k,k);
            A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-A(k+1:n,k)*A(k,k+1:n);
        else
            stop
        end
    end
    %% 分离L=P(LrPr...L1P1)^-1,U
    % P=Pr...P1,Q=Q1...Qr,其中P,Q用不到。
    L = eye(n)+tril(A,-1);
    U = triu(A);
    P=getP_from_u(u,n);
    Q=getQ_from_v(v,n);
    %% 解方程
    %Ax=b , P*A*Q*Q^-1x=P*b , L*U*Q^-1=P*b
    temp=size(u);
    %Ly=Pb
    for i=1:temp(2)
        b=row_trans(i,u(i),b);
    end
    y=SolveLowTriangle(L,b);
    %Uz=y
    z=SolveUpTriangle(U,y);
    %x=Qz
    temp=size(v);
    for i=temp(2):-1:1
        z=row_trans(i,v(i),z);
    end
    x=z;
    % P*A1*Q  =  L*U




    
        
