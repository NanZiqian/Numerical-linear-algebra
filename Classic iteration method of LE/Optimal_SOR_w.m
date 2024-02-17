
% 步长为0.01，0-2，迭代相同次数，取残差更小
function [w] = Optimal_SOR_w(A,b)
    min_error=Inf;
    for omega = 0:0.01:2
        %% 预处理
        N=size(A);
        n=N(1);
        D=diag(A);%D is a column vector containing diag(A)
        L=-tril(A,-1);%L as defined in the textbook
        U=-triu(A,1);%U as defined in the textbook
        M=(diag(D)-omega*L)^-1*( (1-omega)*diag(D)+omega*U );%M is the iterative matrix (D-wL)^-1*[(1-w)D+wU]
        g=left_diag(D.^-1,b);
        B=-left_diag(D.^-1,A);
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
    
    %% 迭代1000
        while 1
            %% 一次迭代 xk=(D-L)^-1U*xk_1+(D-L)^-1*b
            k=k+1;
            xk=repelem(0,n)';% 0 column vector of n dimension
    
            % compute xk
            for i=1:n % compute xk(i)
    
                %求和 \Sigma_{j=1}^{i-1} B(i,j)*xk(j)
                for j=1:i-1 
                    xk(i)=xk(i)+B(i,j)*xk(j);
                end
                %end sum
    
                %求和 \Sigma_{j=i+1}^{n} B(i,j)*xk_1(j)
                for j=i+1:n 
                    xk(i)=xk(i)+B(i,j)*xk_1(j);
                end
                %end sum
                
                xk(i)=omega*(xk(i)+g(i));
                xk(i)=xk(i)+(1-omega)*xk_1(i);
            end
            % end cumpute xk
            if k==10^3
                break
            end
            xk_1=xk;
            
        end % end while
        
        if flag == 0
                error=q/(1-q)*norm(xk-xk_1,inf);
            elseif flag == 1
                error=q/(1-q)*norm(xk-xk_1,1);
            elseif flag == 2
                error=q/(1-q)*norm(xk-xk_1,2); 
        end

        if(error < min_error)
            min_error = error;
            w = omega;
        end
    end %end omega
end

