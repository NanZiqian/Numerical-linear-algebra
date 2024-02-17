
% Jacobi迭代法解线性方程组
function [x] = Jacobi_iteration(A,b)
%% 预处理
    N=size(A);
    n=N(1);
    D=diag(A);%D is a column vector containing diag(A)
    M=left_diag(D.^-1,diag(D)-A);%M is the iterative matrix D^-1(D-A)
    xk_1=norm(b,"inf")*rand(n,1); %x0 initial vector
    k=0;

    q=norm(M,inf);% q is ||M||<1
    flag = 0; % record which norm: 2-norm,2 1-norm,1 inf-norm,0
    if q == 1
        q=norm(M,1);
        flag = 1;
    end
    if q == 1
        q=norm(M,2);
        flag = 2;
    end

    error = 1;% initialize error

%% 迭代至误差范围内
    while error > 10^-8
        %% 一次迭代 xk=D^-1(D-A)*xk_1+D^-1b
        k=k+1;
        xk=repelem(0,n)';% 0 column vector of n dimension

        % compute xk
        for i=1:n % compute xk(i)

            %求和 \Sigma_{j!=i} A(i,j)*xk_1(j)
            for j=1:n 
                if j==i
                    continue
                end
                xk(i)=xk(i)+A(i,j)*xk_1(j);
            end
            %end sum
            
            xk(i)=xk(i)-b(i);
            xk(i)=-1/A(i,i)*xk(i);
        end
        % end cumpute xk
        
        if flag == 0
            error=norm(xk-xk_1,inf);
        elseif flag == 1
            error=norm(xk-xk_1,1);
        elseif flag == 2
            error=norm(xk-xk_1,2); 
        end
        xk_1=xk;
        
    end % end while
    x=xk;

end % end function

