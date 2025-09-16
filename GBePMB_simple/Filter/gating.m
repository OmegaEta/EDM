function [gating_matrix_d,gating_matrix_u,W,ng_in,ng_out] = gating(ppp,mbm,Z,measmodel,paras)
% 基于PPP和MBM的门限，去除杂波。

%ellipsoidal gating for detected targets
%use a boolean vector to store the gating result of each local
%hypothesis
gating_matrix_d = cellfun(@(x1) cell2mat(arrayfun(@(x2) ...
    ellips_gating(x2,Z,measmodel,paras.gating.size),...
    x1,'uniformoutput',false)),mbm.track,'uniformoutput',false);

if paras.ppp_birth
    %ellipsoidal gating for undetected targets
    %use a boolean vector to store the gating result of each component
    gating_matrix_u = cell2mat(arrayfun(@(x) ...
        ellips_gating(x,Z,measmodel,paras.gating.size),...
        ppp,'uniformoutput',false)');
end

%remove unused measurements according to the gating result
%gating result of all targets
%number of tracks
n_track = length(mbm.track);
if paras.mb_birth %whether to use multi-Bernoulli birth model (the default is ppp birth)
    gating = logical(sum(cell2mat(gating_matrix_d'),2));
else
    if n_track > 0
        gating = logical(sum(gating_matrix_u,2) + ...
            sum(cell2mat(gating_matrix_d'),2));
    else
        gating = logical(sum(gating_matrix_u,2));
    end
end

%used measurements
W = Z(:,gating);
%number of measurements after gating
nm = size(W,2);

% 门限内量测数和门限外量测数
ng_in = nm;
ng_out = size(Z,2) - ng_in;

%reconstruct gating matrix
gating_matrix_d = ...
    cellfun(@(x) x(gating,:),gating_matrix_d,'uniformoutput',false);
if ~paras.mb_birth
    gating_matrix_u = gating_matrix_u(gating,:);
end

%% 计算门限内面积
V_ppp = [];
for i = 1:size(ppp,1)
    v = GIWGatingVolume(ppp(i).Cr,ppp(i).v,ppp(i).V,paras,measmodel);
    V_ppp = [V_ppp v];
end
V_mbm = [];
for i = i:size(mbm.track,1)
    v = GIWGatingVolume(mbm.track{i}.Cr,mbm.track{i}.v,mbm.track{i}.V,paras,measmodel);
    V_mbm = [V_mbm v];
end


end

