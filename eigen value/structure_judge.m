
%% 确定最小的l和最大的m，其中H11_l*l为最小的可约上hessenberg矩阵
% H33_m*m为shur标准型
function [l,m]=structure_judge(H)
    %初始化
    [~,n]=size(H);
    l=-1;

    % n <= 2
    if n <= 2
        l=0;
        m=n;
        return
    end

    % n > 2
    % 有连着的次对角线元素不为0
    % case1: m=0
    if H(n,n-1)~=0 && H(n-1,n-2)~=0
        m=0;
    else
    % case2: 1<= m< n-2
        m=1;
        while m < n-2
            if H(n-m+1,n-m) == 0 && H(n-m,n-m-1) ~= 0 && H(n-m-1,n-m-2) ~= 0
                break
            else
                m=m+1;
            end
        end
    % case3: m == n-2
        if m == n-2
            m=n;
            l=0;
        end
    end% end m judging 

    % 寻找到第一个次对角线(l+1,l)为0的l，H11=H(1:l,1:l)
    if l == -1
        for l = n-m:-1:2
            if H(l,l-1) == 0
                l=l-1;
                break
            end
            if l == 2
                l=0;
            end
        end% end for
    end

end
