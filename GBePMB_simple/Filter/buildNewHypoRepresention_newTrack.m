function tracks_new = buildNewHypoRepresention_newTrack(lik_new,bern_new,W,clusters,motionmodel,paras)
%
nm = size(W,2);
dxr = motionmodel.dxr;
n_clusters = size(clusters,2);
if (paras.new_track_init == 2) && paras.mb_approx || paras.mb_birth
    tracks_new = cell(n_clusters,1);
    for c = 1:n_clusters
        tracks_new{c}.c = [0 c];
        tracks_new{c}.lik = [0 lik_new(c)];
        tracks_new{c}.bern = struct('r',0,'xr',zeros(dxr,1),'Cr',...
            ones(dxr,dxr),'V',zeros(2,2),'v',0,'alpha',1,'beta',1,"isDetected",false,'s',paras.Pd_hat*10,'t',(1-paras.Pd_hat)*10,'label',-1,'labellast',-1);
        tracks_new{c}.bern = [tracks_new{c}.bern bern_new(c)];
    end
else
    %create updated single target hypotheses for the first detection of
    %undetected targets
    tracks_new = cell(nm,1);
    for i = 1:nm
        %create single target hypothesis for non-existent target
        tracks_new{i}.c = 0;
        tracks_new{i}.lik = 0;
        tracks_new{i}.bern = struct('r',0,'xr',zeros(dxr,1),'Cr',...
            ones(dxr,dxr),'V',zeros(2,2),'v',0,'alpha',1,'beta',1,"isDetected",false,'s',paras.Pd_hat*10,'t',(1-paras.Pd_hat)*10,'label',-1,'labellast',-1);
    end
    for c = 1:n_clusters
        %create new single target hypothesis
        idx = find(clusters(:,c),1,'last');
        tracks_new{idx}.c = [tracks_new{idx}.c c];% 相当于论文中提出的局部假设表示法，0表示空集
        tracks_new{idx}.lik = [tracks_new{idx}.lik lik_new(c)];
        tracks_new{idx}.bern = [tracks_new{idx}.bern bern_new(c)];
    end
end
end
