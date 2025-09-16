function [AssoTable,IDI] = InitAssociation(priorBern,cell_meas,model)
% 每个KS集关联到每个Bern的似然

N_B = size(priorBern,1);
N_M = size(cell_meas,1);
AssoTable = inf(N_B,N_M);

for i = 1:N_M
    W0 = cell_meas{i}';
    if isempty(W0)
        AssoTable(:,i) = zeros(N_B,1);
        continue;
    end
    for j = 1:N_B
        if priorBern(j).r==0
            AssoTable(j,i) = -inf;
            continue;
        end
        [~,lik] = updateGGIW(priorBern(j).GGIW(end),W0,model);
        N_m = size(W0,2);
        lik = lik/N_m;
        AssoTable(j,i) = lik;
    end
end

IDI = zeros(N_B,N_M);
for i = 1:N_B
   [~,idi] = max(AssoTable(i,:)); 
   IDI(i,idi) = 1;
end%红   绿  蓝 青绿 洋红 黄  黑

end

