
%% 将h_i i-1 <= (|hii|+|h_i-1 i-1|)u 置零
function H=eliminate_Msmall_numbers(H)
    [~,n]=size(H);
    for k = 2:n
        if abs(H(k,k-1)) <= ( abs(H(k,k)) + abs(H(k-1,k-1)) ) * 2^-23
            H(k,k-1)=0;
        end
    end
end