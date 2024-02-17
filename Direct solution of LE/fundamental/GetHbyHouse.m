
%直接从x获得H
function [H] = GetHbyHouse(x)
    [v,beta]=House(x);
    H=GetH_from_v_beta(v,beta);
end

