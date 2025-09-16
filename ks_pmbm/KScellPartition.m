function [Bern_ks_together1,meas_ks_together1] = KScellPartition(N_Bern_prior,N_W,gating_matrix_d1)
% 以gating来对目标和量测的初步分类，以便距离近的目标和量测用来ks划分

Bern_ks_together = cell(N_Bern_prior,1); %可以一起处理量测的bern
meas_ks_together = zeros(N_W,N_Bern_prior);
% Bern_dealed = zeros(N_Bern_prior,1);

gating_matrix_d2 = sum(gating_matrix_d1,2);
gating_matrix_d3 = zeros(N_W,N_Bern_prior);
for i = 1:N_W
    idi = gating_matrix_d1(i,:)>0;
    gating_matrix_d3(i,idi) = gating_matrix_d2(i);
end
for i = 1:N_Bern_prior
    if ~any(gating_matrix_d3(:,i)>1) %如果不同track的门限没有重叠量测
        Bern_ks_together{i} = i;
%         Bern_dealed(i) = 1;
    end
end

%寻找重叠量测集
for i = 1:N_W
   ks_together = find(gating_matrix_d3(i,:)>1);
   for j = ks_together
       Bern_ks_together{j} = [Bern_ks_together{j} ks_together];
   end
end
for i = 1:N_Bern_prior
   Bern_ks_together{i} = unique(Bern_ks_together{i}); 
end
for i = 1:N_Bern_prior
    temp = Bern_ks_together{i};
    for j = Bern_ks_together{i}
         Bern_ks_together{j} = [Bern_ks_together{j} temp]; 
         Bern_ks_together{j} = unique(Bern_ks_together{j});
    end
end
k = 0;
for i = 1:N_Bern_prior
    if ~isempty(Bern_ks_together{i})
        temp = Bern_ks_together{i};
        for j = temp
            meas_ks_together(:,i) = meas_ks_together(:,i) | gating_matrix_d1(:,j);
            if j~=i
                Bern_ks_together{j} = [];
            end
        end
        k = k+1;
        Bern_ks_together1{k,1} = Bern_ks_together{i};
    end
    
end
N_ks_together1 = size(Bern_ks_together1,1);
meas_ks_together1 = cell(N_ks_together1,1);
for i = 1:N_ks_together1
    idi = Bern_ks_together1{i};
    meas_ks_together1{i} = logical(sum(meas_ks_together(:,idi),2));
end