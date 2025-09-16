function [tracks] = falseMeas(W,meas_track,tracks,model)
% 计算量测覆盖的面积，剔除过大的，（这个函数有点简陋）

N_track = size(tracks,1);
N_W = size(W,1);
meas_track1 = cell(N_track,1);
for i = 1:N_track
    N_temp = size(meas_track{i},1);
    meas_track1{i} = false(1,N_W);
    j0 = 0;
    for j = 1:N_temp
        if isempty(meas_track{i}{j})
            continue;
        end
        j0 = j0+1;
        area = GetArea(tracks{i}(j0).Bern.GGIW(end),model);
        idi = find(sum(meas_track{i}{j},1));
        P = area*model.lambda_fa;
        if isempty(idi)
            if P>model.falseMeasT
               N = 0; 
            else
                continue;
            end
        else
            N = size(idi,2);
        end
        
        pp = N/area;
        q = P/pp;
        if q>model.falseMeasT
            tracks{i}(j0).Bern.r = 0;
        end
    end
end


end

