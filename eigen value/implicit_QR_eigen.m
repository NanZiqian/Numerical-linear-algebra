
% 对于任意矩阵，包括拥有复特征值的非对称矩阵；Q记录特征向量，H记录特征值。
function [Q,H] = implicit_QR_eigen(A)
    % 准备工作
    [~,n]=size(A);
    [Q,H]=upper_hessenberg(A);

    while 1
        %收敛性判定
        % 将h_i i-1 <= (|hii|+|h_i-1 i-1|)u 置零
        H=eliminate_Msmall_numbers(H);
        % 确定最小的l和最大的m
        [l,m]=structure_judge(H);
        % 若 m == n，则输出Q,H，分别记录特征向量和特征值
        if m == n
            return
        end

        % double displacement iteration to H22=H(l+1:n-m,l+1:n-m)
        [P,H(l+1:n-m,l+1:n-m)]=double_shift_QR(H(l+1:n-m,l+1:n-m));
        %H=Qk'*H*Qk,H12=H(1:l,l+1:n-m),H23=H(1+l:n-m,n-m+1:n)
        H(1:l,l+1:n-m)=H(1:l,l+1:n-m)*P;
        H(1+l:n-m,n-m+1:n)=P'*H(1+l:n-m,n-m+1:n);
        %Q=Q*Qk
        Q=Q*blkdiag(eye(l),P,eye(m));

    end% end while

end% end function




