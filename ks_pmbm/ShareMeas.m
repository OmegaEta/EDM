function [meas_cell_select,meas_cell_share] = ShareMeas(m,n_tt,ini_likelihood,i_select,meas_cell)
%找出选中量测子集中哪些量测可以被分享。
% 前提：所选量测来自一个目标集
meas_cell_select = meas_cell{i_select};
meas_cell_share = [];
for j = 1:m
    logidi = meas_cell_select==j;
    if ~any(logidi)
        continue;
    end
    [~,idx] = max(ini_likelihood(:,j));
    if i_select==idx && idx<=n_tt && i_select<=n_tt 
        %所选集=置信集&&置信集是一个目标集&&所选量测来自一个目标集
        
    elseif i_select<=n_tt
        %前提：所选量测来自一个目标集
        %不在置信区间的量测可以被分享
        meas_cell_select = meas_cell_select(~logidi);
        meas_cell_share = [meas_cell_share j];
        
    end

end
end


