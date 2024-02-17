

function [x] = SolveLEbyGaussinColumn_max(A,b)
    n = size(A,1);
    %如果输入行向量，则自动转为列向量    
    N=size(b);
    if N(2) ~= 1
        b=b';
    end
%% PA=LU
    for k=1 : n-1
        p=k;
        for i=k+1 : n
            if abs(A(i,k))>abs(A(p,k))
                p=i;
            end
        end
        %交换k行和p行
        A=row_trans(k,p,A);
        v(k)=p;
        if A(k,k) ~= 0
            A(k+1:n,k)=A(k+1:n,k)/A(k,k);
            A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-A(k+1:n,k)*A(k,k+1:n);
        else
            stop
        end
    end
    %% 分离L=P(Ln-1Pn-1...L1P1)^-1,U,P=Pn-1...P1
    L = eye(n)+tril(A,-1);
    U = triu(A);
    P=getP_from_u(v,n);
    %% 解方程
    %Ax=b , PAx=Pb , LUA=Pb
    temp=size(v);
    %Ly=Pb
    for i=1:temp(2)
        b=row_trans(i,v(i),b);
    end
    y=SolveLowTriangle(L,b);
    %Lx=y
    x=SolveUpTriangle(U,y);




    
        
