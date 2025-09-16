function [nc_hat] = gatingBern(mbm,Z,measmodel,paras)
% 
Tc = 1;
[maxw,maxwid] = max(mbm.w);
a_indices = mbm.table(maxwid,:);
n_track = length(mbm.track);
if n_track>0
    single_hypo = [];
    for j = 1:n_track
        if a_indices(j) > 0
            single_hypo = [single_hypo;mbm.track{j}(a_indices(j))];
        end
    end
    
    gating_matrix_d = arrayfun(@(x) ellips_gating(x,Z,measmodel,paras.gating.size),single_hypo,'uniformoutput',false);
    gating_matrix = sum(cell2mat(gating_matrix_d'),2)>0;
    n = size(Z,2);
    ng_in = sum(gating_matrix,1);
    nc_out_hat = n - ng_in;
    
    % 
    v = cell2mat(arrayfun(@(x) GIWGatingVolume(x.Cr,x.v,x.V,paras,measmodel),single_hypo,'uniformoutput',false));
    Vs = sum(v,1);
    psg = measmodel.c_pdf * Vs;
    
    % 目标量测率
    Lambda_X = cell2mat(arrayfun(@(x) (x.r>0.5)*x.isDetected*x.alpha/x.beta,single_hypo,'uniformoutput',false));
    lambda_X = sum(Lambda_X,1);

    w_out = 1 - psg;
    w_in = psg;
    nc_in_hat = optimizeClutter(ng_in,psg,lambda_X,Tc);
    nc_hat = nc_out_hat*w_out+nc_in_hat*w_in;
else
    nc_hat = size(Z,2);
end





end

