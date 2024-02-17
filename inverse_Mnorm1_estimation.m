function [v] = inverse_Mnorm1_estimation(A)
%INVERSE_MNORM1_ESTIMATION 此处显示有关此函数的摘要
%   此处显示详细说明
    n=size(A,1);
    for k=1:n
        x(k,1)=1/k;
    end
    while 1
        x1=x;
        x1=SolveLEbyGaussinColumn_max(A',x);
        sgn=sign(x1);
        sgn=SolveLEbyGaussinColumn_max(A,sgn);
        if norm(sgn,"inf")<sgn'*x
            v=norm(x1,1);
            break;
        else 
            i=findmax_v(sgn);
            x=repelem(0,n);
            x(i)=1;
            x=x';
        end
    end
end

