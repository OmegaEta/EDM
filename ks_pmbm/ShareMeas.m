function [meas_cell_select,meas_cell_share] = ShareMeas(m,n_tt,ini_likelihood,i_select,meas_cell)
%�ҳ�ѡ�������Ӽ�����Щ������Ա�����
% ǰ�᣺��ѡ��������һ��Ŀ�꼯
meas_cell_select = meas_cell{i_select};
meas_cell_share = [];
for j = 1:m
    logidi = meas_cell_select==j;
    if ~any(logidi)
        continue;
    end
    [~,idx] = max(ini_likelihood(:,j));
    if i_select==idx && idx<=n_tt && i_select<=n_tt 
        %��ѡ��=���ż�&&���ż���һ��Ŀ�꼯&&��ѡ��������һ��Ŀ�꼯
        
    elseif i_select<=n_tt
        %ǰ�᣺��ѡ��������һ��Ŀ�꼯
        %�������������������Ա�����
        meas_cell_select = meas_cell_select(~logidi);
        meas_cell_share = [meas_cell_share j];
        
    end

end
end


