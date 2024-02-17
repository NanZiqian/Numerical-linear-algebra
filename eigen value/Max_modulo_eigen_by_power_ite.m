
%  幂法求模最大特征值与其特征向量
function [lambda1,v1,flag,lambda2,v2] = Max_modulo_eigen_by_power_ite(A)
%MAX_MODULO_EIGEN 此处显示有关此函数的摘要
    n=size(A,1);
    uk=zeros(n,3);
    muk=[0,0];

    % 初始化uk,uk(:,3),uk(:,2),uk(:,1)分别代表uk_2,uk_1,uk
    uk(:,3)=rand(n,1);
    uk(:,3)=uk(:,3)/norm(uk(:,3),"inf");
    it=0;

    % 迭代主体,muk(2)代表muk,muk(1)代表muk_1
    while norm(uk(:,3)-uk(:,2),"inf")>10^-8 && norm(uk(:,3)-uk(:,1),"inf")>10^-8
        it=it+1;
        uk(:,1)=[];% uk_1=uk,uk_2=uk_1;
        muk(1)=[];% muk_1=muk

        yk=A*uk(:,2);% yk=A*uk_1;

        i=find_max_modulo_index(yk);
        muk(2)=yk(i);%yk模最大的分量
        %muk(2) = norm(yk,"inf");

        uk(:,3)=yk/muk(2);% uk=yk/muk;

        %break condition
        %如果需要更多的判断条件，while 1 ，在这里if break 就行
        if it > 10^5
            disp('Max_modulo_eigen_by_power_ite: it reaches 10^5!');
            break
        end
    end

    %% 判断矩阵特征值分布
    % flag = 1: |lambda1| > |lambda2|
    % flag = 2: lambda1 = -lambda2
        % case 1: |lambda1| > |lambda2| , lambda1 > 0
    if norm(uk(:,3)-uk(:,2),"inf") <= 10^-8
        v1 = uk(:,3);
        lambda1 = muk(2);
        flag = 1;
        lambda2 = Inf;
        v2 = repelem(Inf,n)';
    elseif abs(uk(:,3)+uk(:,2)) < 10^-7
        % case 2: |lambda1| > |lambda2| , lambda1 < 0
        %如果用模最大分量，则lambda1 < 0 uk奇偶数列也共极限，这里用不到
        v1 = uk(:,3);
        lambda1 = -muk(2);
        flag = 1;
        lambda2 = Inf;
        v2 = repelem(Inf,n)';
    elseif norm(uk(:,3)-uk(:,1),"inf") <= 10^-8
        % case 3: lambda1 = -lambda2
        lambda1 = sqrt(muk(2)*muk(1));
        v1 = uk(:,2)*sqrt(muk(1)/muk(2))+uk(:,1);
        flag = 2;

        lambda2 = -lambda1;
        v2 = -uk(:,3)*sqrt(muk(2)/muk(1))+uk(:,2);
    end% 如果要更多的条件，可以继续加elseif
end

