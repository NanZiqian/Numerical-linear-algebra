
%这是从v和b得到Householder矩阵H的函数
function [H] = GetH_from_v_beta(v,beta)
    %如果输入列向量，则自动转为行向量    
    N=size(v);
    if N(1) ~= 1
        v=v';
        n=N(1);
    else
        n=N(2);
    end
    %获得H
    H=eye(n)-beta*v'*v;
end

