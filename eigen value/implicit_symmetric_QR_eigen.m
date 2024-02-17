
%   此处显示详细说明
function [Q,T] = implicit_symmetric_QR_eigen(A)
    [~,n]=size(A);
    [Q,T]=tri_diagonalize(A);

    while 1
        %收敛性判定
        % 将h_i i-1 <= (|hii|+|h_i-1 i-1|)u 置零
        T=eliminate_Msmall_numbers(T);
        T=T';
        T=eliminate_Msmall_numbers(T);
        T=T';
        % 确定最小的l和最大的m
        [l,m]=structure_judge(T);
        % 若 m == n，则输出Q,H，分别记录特征向量和特征值
        if m == n
            return
        end

        [P,T(l+1:n-m,l+1:n-m)]=wilkinson_shift_symmetric_QR(T(l+1:n-m,l+1:n-m));
        Q=Q*blkdiag(eye(l),P,eye(m));
    end

end

