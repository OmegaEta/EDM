function [lik_new,bern_new] = everyClusterIsNewBern(ppp,W,clusters,gating_matrix_u,motionmodel,measmodel,paras)
%
dxr = motionmodel.dxr;
n_clusters = size(clusters,2);
bern_new = repmat(struct('r',0,'xr',zeros(dxr,1),'Cr',ones(dxr,dxr),...
    'V',zeros(2,2),'v',0,'alpha',1,'beta',1,'isDetected',false,'s',paras.Pd_hat*10,'t',(1-paras.Pd_hat)*10,'label',-1,'labellast',-1),1,n_clusters);
lik_new = zeros(n_clusters,1);
if paras.mb_birth
    for c = 1:n_clusters
        lik_new(c) = sum(clusters(:,c))*log(paras.clutter.c_intensity_hat);
    end
else
    for c = 1:n_clusters
        %check if the cth cluster is in the gate of any ppp components
        ppp_idx = sum(gating_matrix_u-clusters(:,c)<0) == 0;
        if any(ppp_idx)
            [lik_new(c),bern_new(c)] = ...
                BGGppp_upd(ppp(ppp_idx),W(:,clusters(:,c)),measmodel,paras);
        else
            %if not, they are all clutter
            lik_new(c) = sum(clusters(:,c))*log(paras.clutter.c_intensity_hat);
        end
    end
end
end