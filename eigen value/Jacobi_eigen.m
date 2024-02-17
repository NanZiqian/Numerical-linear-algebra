
%   经典Jacobi迭代,Qk^H*A*Qk=Ak，其中Ak非对角元很小
function [Q,A] = Jacobi_eigen(A)
    if A ~= A'
        dips('Jacobi_eigen: A is not symmetric!')
        return 
    end
    [~,n]=size(A);
    Q=eye(n);
    delta_sqaure = (10^-8*norm(A,'fro'))^2;
    E_A_square=norm(A,'fro')^2-norm(diag(A),2)^2;
    while E_A_square > delta_sqaure
        [p,q]=max_module_not_diagonal_element(A);
        [c,s,t]=GivensInJacobi(A,p,q);

        %% A = G(p,q,theta)^H*A*G(p,q,theta)
        for i=1:n
            if i == p || i == q
                continue
            end

            % A(i,p)=A(p,i)=beta_ip(i)
            beta_ip(i)=c*A(i,p)-s*A(i,q);
            % A(i,q)=A(q,i)=beta_iq(i)
            beta_iq(i)=s*A(i,p)+c*A(i,q);
        end
        
        beta_pp=c^2*A(p,p)-2*s*c*A(p,q)+s^2*A(q,q);
        beta_qq=s^2*A(p,p)+2*s*c*A(p,q)+c^2*A(q,q);

        %check correctness
        beta_pq=(c^2-s^2)*A(p,q)+s*c*(A(p,p)-A(q,q));
        if abs(beta_pq) > 10^-8
            disp('Jacobi_eigen: beta_pq is not 0!');
            return
        end

        for i=1:n
            if i == p || i == q
                continue
            end
            A(p,i)=beta_ip(i);
            A(i,p)=A(p,i);
            A(q,i)=beta_iq(i);
            A(i,q)=A(q,i);
        end
        A(p,p)=beta_pp;
        A(q,q)=beta_qq;
        A(p,q)=0;
        A(q,p)=0;
        % end A = G(p,q,theta)^H*A*G(p,q,theta)
        
        G=eye(n);
        G(p,p)=c;
        G(q,q)=c;
        G(p,q)=s;
        G(q,p)=-s;
        Q=Q*G;

        E_A_square=norm(A,'fro')^2-norm(diag(A),2)^2;
    end % end while
end

function [c,s,t]=GivensInJacobi(A,p,q)
    if A(p,q) ~= 0
        tau = (A(q,q)-A(p,p))/2/A(p,q);
        if tau >= 0
            t = 1/(tau+ sqrt( 1+tau^2 ) );
        else
            t = 1/( tau - sqrt( 1+tau^2 ) );
        end
        c = 1/sqrt( 1+t^2 );
        s = t*c;
    else
        c=1;s=0;t=0;
    end
end