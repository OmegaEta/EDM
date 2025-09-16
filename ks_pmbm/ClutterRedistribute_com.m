function [cell_meas,C_t] = ClutterRedistribute_com(cell_splited,N,C_f,model)
% clarify
d = 2;
P_C = 0;
N_C = size(C_f,1);
if N==0
    C_t = [];
    cell_meas{1}=C_f;
    return;
end
adaGauss = cell(N,1);
for i = 1:N
    [adaGauss{i}] = adaMeas(cell_splited{i},model);
end
selected = false(N_C,1);
%%
%Ëæ»ú²ÉÑù
for i = 1:N_C
    i_selected = i;
    selected(i_selected) = true;
    PDF = [];
    for j = 1:N
        Sigma = adaGauss{j}.P;
        mu = adaGauss{j}.m;
        pdf = mvnpdf(C_f,mu,Sigma)';
        PDF = [PDF;pdf];
    end

    P_N = sum(PDF,1);
    P_k = PDF./repmat(P_N,N,1);
    [~,i_combine] = max(P_k(:,i_selected));
    cell_splited{i_combine} = [cell_splited{i_combine};C_f(i_selected,:)];
    adaGauss = cell(N,1);
    for i0 = 1:N
        [adaGauss{i0}] = adaMeas(cell_splited{i0},model);
    end
end

%%
C_t = cell_splited{N+1};
cell_meas = cell(N,1);
for i = 1:N
    cell_meas{i} = cell_splited{i};
end
end

