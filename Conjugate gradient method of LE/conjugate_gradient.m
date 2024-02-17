
% CONJUGATE_GRADIENT，共轭梯度法求对称正定方程组
% x_0为初值，可任意非0,此处取rand(n,1)
function [x] = conjugate_gradient(A,b)
    
    N=size(A);
    n=N(1);
    k=0;
    % xk_1=rand(n,1);
    xk_1=[0;0;1]
    rk_1=b-A*xk_1;
    rk=1;
    while norm(abs(rk),1)>10^-8
        k = k+1;
        %% begin kth iteration,after we'll get xk
        if k == 1
            %calculate p_{k-1}, beta_{k-1}=0
            pk_1=rk_1;
        else
            
            %update coefficient xk_1,rk_2,rk_1,pk_2
            xk_1=xk;

            rk_2=rk_1;
            rk_1=rk;

            pk_2=pk_1;
            %

            %calculate beta_{k-1} and p_{k-1}
            betak_1=rk_1'*rk_1/(rk_2'*rk_2);
            pk_1=rk_1+betak_1*pk_2;

        end
        alpha=rk_1'*rk_1/(pk_1'*A*pk_1);
        xk=xk_1+alpha*pk_1;
        rk=rk_1-alpha*A*pk_1;
    end % end while
    x=xk;
    
end

