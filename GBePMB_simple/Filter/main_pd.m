clear;
rng('default');

addpath('..\clutterEstimate\');
addpath('..\scenarioCreate\');
addpath("..\MEOT-GGIW\");
addpath('..\Third-Party-Code\');
addpath('E:\matlab files\研究生\Coexisting-point-extended-target-PMBM-filter-main\Coexisting-point-extended-target-PMBM-filter-main\Third-party code');

%% 模型参数初始化
motionmodel = paraOfMotionmodel(1);
birthmodel = paraOfBirthmodel(4);
measmodel = paraOfMeasmodel(1);
%% 产生真实轨迹
T=100;
[gt,gt_time,card] = generateGroundtruth(2,birthmodel,motionmodel,T);
nt = size(gt,2);
%% 产生量测
Z = generateMeas(gt,measmodel,motionmodel,T);
%% 绘制真实轨迹，绘制量测
plot_enable_trajectory = 0;
plot_enable_result = 1;
if plot_enable_trajectory
    pt = 0.01;
    plotTrajectory(gt,motionmodel,1);
    plotDynamicTrajectory(gt_time,gt,motionmodel,2,pt);
    % plotDynamicMeas(Z,3,pt);
end
%% 滤波器参数初始化
[paras,ppp,mbm] = parasSetting(birthmodel,measmodel,2);

%memory to store state estimate
est = cell(T,1);
card_est = zeros(T,1);
t_elapsed = zeros(T,1);

%parameters of GOSPA metric
gospa.p = 1;
gospa.c = 20;
gospa.alpha = 2;

lambdaC = [];
Nc_hat = [];
%recursive Bayesian estimation
fprintf('Time step: ')

global Label_id;
Label_id = 0;
%% 滤波
if plot_enable_trajectory
    figure(8);
    clf(8);
    grid on
    axis([-200,200,-200,200]);
    hold on;
end
for i = 1:nt
    %plot trajectory
    gt_ks_plot = plot(gt(i).xr(1,:),gt(i).xr(2,:),'b','linewidth',0.3);
end
for t = 1:T
    fprintf('%d ',t)
    tic
    if plot_enable_trajectory
        H = [];
        h_Z = plot(Z{t}(1,:),Z{t}(2,:),"r*");
        H = [H h_Z];
    end
    % 门限
    [gating_matrix_d,gating_matrix_u,W,ng_in,ng_out] = gating(ppp,mbm,Z{t},measmodel,paras);
    nm = size(W,2);
    n_track = length(mbm.track);
    
    %% 划分
    [partitions,clusters,np,n_clusters,partitions_indices] = partition(W,mbm,paras);
    %% 假设每个簇都是一个新的Bern
    [lik_new,bern_new] = everyClusterIsNewBern(ppp,W,clusters,gating_matrix_u,motionmodel,measmodel,paras);

    %% 存在的track -- 簇 更新
    %create updated single target hypotheses for detected targets
    %number of single target hypotheses per track
    tracks_upd = existTracksUpdate(mbm,W,clusters,gating_matrix_d,measmodel,paras);
 
    %% 构建O(n)级别的假设表示法（n为量测数）
    tracks_new = buildNewHypoRepresention_newTrack(lik_new,bern_new,W,clusters,motionmodel,paras);
    
    %% 把上个时刻的mbm_upd，更新为本时刻的mbm_upd
    %reset m to the number of new tracks
    m = length(tracks_new);
    
    %update local hypothesis trees
    mbm_upd.track = cell(n_track+m,1);
    for i = 1:n_track
        idx = 0;
        for j = 1:length(tracks_upd{i})
            mbm_upd.track{i} = [mbm_upd.track{i} tracks_upd{i}{j}.bern];
            %use an extra variable to record the index of each new single
            %target hypothesis in local hypothesis tree i
            n_ij = length(tracks_upd{i}{j}.c);
            tracks_upd{i}{j}.idx = (1:n_ij) + idx;
            idx = idx + n_ij;
        end
    end
    for i = 1:m
        mbm_upd.track{n_track+i,1} = tracks_new{i}.bern;
    end
    %%
    %data association for each global association hypothesis
    mbm_upd.w = [];
    mbm_upd.table = zeros(0,n_track+m);
    A = length(mbm.w);
    % 
    for a = 1:A
        a_indices = mbm.table(a,:);
        costs_temp = 0;
        for j = 1:n_track
            if a_indices(j) > 0
                single_hypo = tracks_upd{j}{a_indices(j)};
                costs_temp = costs_temp + single_hypo.lik(1);%costs_temp为当前全局假设中的所以局部假设都miss
            end
        end
        %if there is no measurement partition, all the targets are
        %misdetected
        if n_clusters==0
            table_upd = zeros(1,n_track+m);
            w_upd = 0;
            for j = 1:n_track
                if a_indices(j) > 0
                    table_upd(1,j) = tracks_upd{j}{a_indices(j)}.idx(1);
                    w_upd = w_upd + tracks_upd{j}{a_indices(j)}.lik(1);
                end
            end
            for j = n_track+1:n_track+m
                table_upd(1,j) = 1;
            end
            mbm_upd.w = [mbm_upd.w;w_upd + mbm.w(a)];
            mbm_upd.table = [mbm_upd.table;table_upd];
        end
        %go through each measurement partition
        for p = 1:np
            %find all the clusters under this partition
            p_idx = partitions_indices{p};
            %number of clusters under this partition
            p_n = length(p_idx);
            %construct cost matrix
            C = inf(p_n,n_track+m);
            for j = 1:n_track
                if a_indices(j) > 0
                    single_hypo = tracks_upd{j}{a_indices(j)};
                    %set cost for detected targets
                    [LIA,LOCB] = ismember(single_hypo.c,p_idx);
                    C(LOCB(LOCB>0),j) = -single_hypo.lik(LIA)' + single_hypo.lik(1);
                end
            end
            %set cost for undetected targets
            for j = n_track+1:n_track+m
                single_hypo = tracks_new{j-n_track};
                [LIA,LOCB] = ismember(single_hypo.c,p_idx);
                C(LOCB(LOCB>0),j) = -single_hypo.lik(LIA)';
            end
            
            %find columns of C that contain finite entries
            idx = find(sum(isfinite(C),1) > 0);
            C = C(:,idx);
            if paras.data_assoc == 1
                %find M-best assignments using Murty's algorithm
                [assignments,~,costs] = kBest2DAssign(C,ceil(paras.M*exp(mbm.w(a))));
                assignments = assignments';
                costs = costs';
            elseif paras.data_assoc == 2
                %find M-best assignments using Gibbs sampling
                [assignments,costs] = assign2DByGibbs(C,...
                    paras.gibbs,paras.M);
            end
            %number of assignments
            n_a = size(assignments,1);
            %restore track indices
            for i = 1:n_a
                assignments(i,:) = idx(assignments(i,:));
            end
            
            table_upd = zeros(n_a,n_track+m);
            %update the glocal hypothesis look-up table and weight
            for i = 1:n_a
                %go through each association in a given assignment
                for j = 1:p_n
                    %check if detected targets or undetected targets
                    assoc = assignments(i,j);
                    if assoc <= n_track
                        %find the index of the corresponding cluster
                        table_upd(i,assoc) = tracks_upd{assoc}...
                            {a_indices(assoc)}.idx(tracks_upd{assoc}...
                            {a_indices(assoc)}.c == p_idx(j));
                    else
                        %find the index of the corresponding cluster
                        table_upd(i,assoc) = ...
                            find(tracks_new{assoc-n_track}.c == p_idx(j),1);
                    end
                end
                %go through each unassociated track
                unassign = true(n_track+m,1);
                unassign(assignments(i,:)) = false;
                unassign = find(unassign);
                for j = 1:length(unassign)
                    if unassign(j) <= n_track
                        temp = a_indices(unassign(j));
                        if temp > 0
                            table_upd(i,unassign(j)) = ...
                                tracks_upd{unassign(j)}{temp}.idx(1);
                        else
                            table_upd(i,unassign(j)) = 0;
                        end
                    else
                        table_upd(i,unassign(j)) = 1;
                    end
                end
            end
            
            %when computing the global hypothesis weight, assume that each
            %measurement partition is equal likely, i.e., uniform
            %distributed
            mbm_upd.w = [mbm_upd.w;-costs'+costs_temp+mbm.w(a)];
            mbm_upd.table = [mbm_upd.table;table_upd];
        end
    end
    
    %normalise global hypothesis weight
    mbm_upd.w = normalizeLogWeights(mbm_upd.w);
    %%
    %prune updated global hypotheses with small weights
    [mbm_upd_w,order] = sort(exp(mbm_upd.w),'descend');
    pos = find(cumsum(mbm_upd_w)>=1-exp(paras.pruning.w),1);
    mbm_upd.w = mbm_upd.w(order(1:pos));
    mbm_upd.table = mbm_upd.table(order(1:pos),:);
    mbm_upd.w = normalizeLogWeights(mbm_upd.w);
    if ~paras.mb_approx
        %cap the number of global hypotheses
        if length(mbm_upd.w) > paras.cap.w
            [~,idx] = sort(mbm_upd.w,'descend');
            mbm_upd.w = mbm_upd.w(idx(1:paras.cap.w));
            mbm_upd.w = normalizeLogWeights(mbm_upd.w);
            mbm_upd.table = mbm_upd.table(idx(1:paras.cap.w),:);
        end
    end
    
    %remove single target hypotheses with small probability of existence
    for i = 1:length(mbm_upd.track)
        if paras.mb_birth
            idx = find(([mbm_upd.track{i}.r] < paras.pruning.r) | ...
                ([mbm_upd.track{i}.alpha]./[mbm_upd.track{i}.beta] < 1) | ...
                ([mbm_upd.track{i}.v] < 7));
        else
            idx = find(([mbm_upd.track{i}.r] < exp(paras.pruning.ppp)) | ...
                ([mbm_upd.track{i}.alpha]./[mbm_upd.track{i}.beta] < 1) | ...
                ([mbm_upd.track{i}.v] < 7));
        end
        for j = 1:length(idx)
            mbm_upd.table(mbm_upd.table(:,i) == idx(j),i) = 0;
        end
        %re-index
        idx_0 = mbm_upd.table(:,i) > 0;
        [idx,~,temp] = unique(mbm_upd.table(idx_0,i));
        mbm_upd.table(idx_0,i) = temp;
        mbm_upd.track{i} = mbm_upd.track{i}(idx);
    end
    %remove empty track
    idx = ~cellfun('isempty',mbm_upd.track);
    mbm_upd.track = mbm_upd.track(idx);
    mbm_upd.table = mbm_upd.table(:,idx);
    if isempty(mbm_upd.table)
        mbm_upd.table = zeros(1,0);
        mbm_upd.track = cell(0,1);
        mbm_upd.w = 0;
    end
    
    %merge rows of global hypothesis look-up table that are the same
    if length(mbm_upd.w) > 1
        [mbm_upd.table,~,IC] = unique(mbm_upd.table,'rows');
        n_a = size(mbm_upd.table,1);
        temp = zeros(n_a,1);
        for i = 1:n_a
            [~,temp(i)] = normalizeLogWeights(mbm_upd.w(IC==i));
        end
        mbm_upd.w = temp;
    end
    
    %misdetection update of ppp
    if ~paras.mb_birth
        ppp = BGGppp_miss(ppp,measmodel,paras);
    end
    
    %number of global hypotheses and tracks
    n_mb = size(mbm_upd.table,1);
    %%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if paras.mb_approx && n_mb > 1
        %find tracks with only one single target hypothesis being included
        %in any of the global hypotheses
        idx = sum(mbm_upd.table - 1,1) ~= 0;
        mbm_upd.track = [mbm_upd.track(idx);mbm_upd.track(~idx)];
        mbm_upd.table = [mbm_upd.table(:,idx) mbm_upd.table(:,~idx)];
        
        %track-oriented merging
        %number of tracks with conflicts
        n_track_t = length(idx);
        for i = 1:n_track_t
            %compute marginal probability that the single target hypothesis
            %in track i is included in the global hypothesis
            c = unique(mbm_upd.table(:,i));
            c = c(c>0);
            lc = length(c);
            w_margin = zeros(lc,1);
            for j = 1:lc
                %find all the global hypotheses that contain this single
                %target hypothesis and compute the log sum of the weights
                [~,log_sum_w] = normalizeLogWeights(mbm_upd.w(mbm_upd.table(:,i)==c(j)));
                %also take into account the probability of existence
                w_margin(j) = log_sum_w + log(mbm_upd.track{i}(c(j)).r);
            end
            %perform moment matching
            [w_margin_n,log_sum_w] = normalizeLogWeights(w_margin);
            mb_approx(i).r = exp(log_sum_w);
            [mb_approx(i).xr,mb_approx(i).Cr] = ...
                kinematic_merge(mbm_upd.track{i},w_margin_n);
            [mb_approx(i).V,mb_approx(i).v] = ...
                extent_merge(mbm_upd.track{i},w_margin_n);
            [mb_approx(i).alpha,mb_approx(i).beta] = ...
                gamma_merge(mbm_upd.track{i},w_margin_n);
            mb_approx(i).isDetected = ...
                detected_merge(mbm_upd.track{i},w_margin_n);
            [mb_approx(i).s,mb_approx(i).t] = ...
                beta_merge1(mbm_upd.track{i},w_margin_n);
            mb_approx(i).label = mbm_upd.track{i}(1).label;
        end
        
        if paras.mb_merge == 1 || n_track_t < 2
            %track-oriented merging
            for i = 1:n_track_t
                mbm_upd.track{i} = mb_approx(i);
            end
        elseif paras.mb_merge == 2 && n_track_t > 1
            %initialized using the result of track-oriented merging
            %only need to consider track with more than one single target
            %hypothesis
            bern0 = struct('r',0,'xr',zeros(dxr,1),...
                'Cr',ones(dxr,dxr),'V',zeros(2,2),'v',0,...
                'alpha',1,'beta',1);
            mbm_hypo = repmat(bern0,[n_mb,n_track_t]);
            for i = 1:n_mb
                for j = 1:n_track_t
                    if mbm_upd.table(i,j) > 0
                        mbm_hypo(i,j) = mbm_upd.track{j}(mbm_upd.table(i,j));
                    end
                end
            end
            nh = cellfun('length',mbm_upd.track(1:n_track_t));
            
            %iterative optimization
            C_optimal = inf;
            mb_hat = mb_approx;
            while(1)
                Cmin = zeros(n_mb,1);
                assignment_mb = zeros(n_track_t,n_mb);
                
                Ch = cell(n_track_t,n_track_t);
                for k = 1:n_track_t
                    for i = 1:n_track_t
                        Cht = zeros(nh(i),1);
                        for j = 1:nh(i)
                            Cht(j) = bern_cross_entropy(mbm_upd.track{i}(j),mb_hat(k));
                        end
                        Ch{k,i} = Cht;
                    end
                end
                
                for i = 1:n_mb
                    %construct cost matrix with columns corresponding
                    %to Bernoullis components in the ith MB and rows
                    %corresponding to the approximated Bernoullis
                    C = zeros(n_track_t,n_track_t);
                    for j = 1:n_track_t
                        for k = 1:n_track_t
                            if mbm_upd.table(i,j) == 0
                                C(k,j) = bern_cross_entropy(bern0,mb_hat(k));
                            else
                                C(k,j) = Ch{k,j}(mbm_upd.table(i,j));
                            end
                        end
                    end
                    %find the optimal assignment, assignment_mb(j,i)
                    %means that the jth approximated Bernoulli should
                    %be merged using the assignment_mb(j,i)th Bernoulli
                    %in the ith MB
                    [assignment_mb(:,i),~,Cmin(i)] = assign2D(C);
                end
                %check the columns of assignment_mb, if they are all of
                %the form 1:n_b, then the iterative optimization
                %converges
                if ~any(assignment_mb - (1:n_track_t)')
                    break;
                else
                    C_hat = Cmin'*exp(mbm_upd.w);
                    if C_optimal-C_hat < paras.vb_threshold
                        break;
                    else
                        C_optimal = C_hat;
                    end
                end
                %if not converged, perform merging to construct new
                %mb_hat
                for j = 1:n_track_t
                    w_margin = zeros(n_mb,1);
                    bern = repmat(bern0,[1,n_mb]);
                    for i = 1:n_mb
                        w_margin(i) = mbm_upd.w(i) + ...
                            log(mbm_hypo(i,assignment_mb(j,i)).r);
                        bern(i) = mbm_hypo(i,assignment_mb(j,i));
                    end
                    [w_margin_n,log_sum_w] = normalizeLogWeights(w_margin);
                    %note that each MB may contain Bernoullis with zero
                    %probability of existence with unvalid pdf. The
                    %merging of these Bernoullis will of course be a
                    %Bernoulli with zero probability of existence
                    if isnan(w_margin_n)
                        mb_hat(j) = bern0;
                    else
                        mb_hat(j).r = exp(log_sum_w);
                        [mb_hat(j).xr,mb_hat(j).Cr] = kinematic_merge(bern,w_margin_n);
                        [mb_hat(j).V,mb_hat(j).v] = extent_merge(bern,w_margin_n);
                        [mb_hat(j).alpha,mb_hat(j).beta] = gamma_merge(bern,w_margin_n);
                        mb_hat(i).isDetected = detected_merge(bern_h,w_margin_n);
                        [mb_hat(j).s,mb_hat(j).t] = beta_merge1(bern,w_margin_n);
                    end
                end
            end
            for i = 1:n_track_t
                mbm_upd.track{i} = mb_hat(i);
            end
            
        elseif paras.mb_merge == 3 && n_track_t > 1
            %initialized using the result of track-oriented merging
            %only need to consider track with more than one single target
            %hypothesis
            bern0 = struct('r',0,'xr',zeros(dxr,1),...
                'Cr',ones(dxr,dxr),'V',zeros(2,2),'v',0,...
                'alpha',1,'beta',1,'isDetected',false,'s',2,'t',2,"label",0);
            nh = cellfun('length',mbm_upd.track(1:n_track_t))+1;
            bern_h = [];
            for i = 1:n_track_t
                bern_h = [bern_h bern0];
                bern_h = [bern_h mbm_upd.track{i}];
            end
            
            n_H = sum(nh);
            %construct constraints
            qhj = zeros(n_H,n_track_t);
            idxi = 0;
            for i = 1:n_track_t
                for j = 1:nh(i)
                    qhj(j+idxi,i) = sum(exp(mbm_upd.w((mbm_upd.table(:,i) == j-1))));
                end
                idxi = idxi + nh(i);
            end
            idx = sum(qhj,2)>0;
            bern_h = bern_h(idx);
            n_H = length(bern_h);
            qhj = qhj(idx,:);
            qj = sum(qhj,1)';
            qh = sum(qhj,2);
            
            %iterative optimization
            C_optimal = inf;
            mb_hat = mb_approx;
            while(1)
                %construct cost matrix
                C = zeros(n_H,n_track_t);
                for i = 1:n_H
                    for j = 1:n_track_t
                        C(i,j) = bern_cross_entropy(bern_h(i),mb_hat(j));
                    end
                end
                %solve the transportation problem
                [C_hat,q_hat] = LP_transport(C,qh,qj);
                %if q_hat and qhj are very similar, then the iterative
                %optimization converges
                if sum(sum(abs(q_hat-qhj))) < 1e-6
                    break;
                end
                if C_optimal-C_hat < paras.vb_threshold
                    break;
                else
                    C_optimal = C_hat;
                end
                %if not converged, perform merging to construct new
                %mb_hat
                for i = 1:n_track_t
                    mb_hat(i).r = [bern_h.r]*q_hat(:,i);
                    if mb_hat(i).r > 0
                        w_margin_n = log([bern_h.r]'.*q_hat(:,i)/mb_hat(i).r);
                        [mb_hat(i).xr,mb_hat(i).Cr] = kinematic_merge(bern_h,w_margin_n);
                        [mb_hat(i).V,mb_hat(i).v] = extent_merge(bern_h,w_margin_n);
                        [mb_hat(i).alpha,mb_hat(i).beta] = gamma_merge(bern_h,w_margin_n);
                        mb_hat(i).isDetected = detected_merge(bern_h,w_margin_n);
                        [mb_hat(j).s,mb_hat(j).t] = beta_merge1(bern,w_margin_n);
                    else
                        mb_hat(i) = bern0;
                    end
                end
            end
            for i = 1:n_track_t
                mbm_upd.track{i} = mb_hat(i);
            end
        end
        mbm_upd.w = 0;
        mbm_upd.table = ones(1,n_track_t);
        
        %recycle tracks with small probability of existence
        idx_logical = (cellfun(@(x) x.r, mbm_upd.track) >= paras.pruning.r) & ...
            (cellfun(@(x) x.alpha, mbm_upd.track)./...
            cellfun(@(x) x.beta, mbm_upd.track) > 1) & ...
            (cellfun(@(x) x.v, mbm_upd.track) > 7);
        if paras.recycle && (~paras.mb_birth)
            idx = find(~idx_logical);
            for j = 1:length(idx)
                ppp(end+1,1).w = log(mbm_upd.track{idx(j)}.r);
                ppp(end,1).xr = mbm_upd.track{idx(j)}.xr;
                ppp(end,1).Cr = mbm_upd.track{idx(j)}.Cr;
                ppp(end,1).V = mbm_upd.track{idx(j)}.V;
                ppp(end,1).v = mbm_upd.track{idx(j)}.v;
                ppp(end,1).alpha = mbm_upd.track{idx(j)}.alpha;
                ppp(end,1).beta = mbm_upd.track{idx(j)}.beta;
                ppp(end,1).isDetected = false;
                ppp(end,1).s = mbm_upd.track{idx(j)}.s;
                ppp(end,1).t = mbm_upd.track{idx(j)}.t;
            end
        end
        %delete these single target hypotheses
        mbm_upd.track = mbm_upd.track(idx_logical);
        mbm_upd.table = mbm_upd.table(idx_logical);
        if isempty(mbm_upd.table)
            mbm_upd.table = zeros(1,0);
            mbm_upd.track = cell(0,1);
            mbm_upd.w = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        %recycle single target hypotheses with small probability of existence
        for i = 1:length(mbm_upd.track)
            idx = find(([mbm_upd.track{i}.r] < paras.pruning.r) | ...
                ([mbm_upd.track{i}.alpha]./[mbm_upd.track{i}.beta] < 1) | ...
                ([mbm_upd.track{i}.v] < 7));
            %for each single target hypothesis to be recycled, find global
            %hypotheses that include it
            for j = 1:length(idx)
                temp = mbm_upd.table(:,i)==idx(j);
                mbm_upd.table(temp,i) = 0;
                if paras.recycle && (~paras.mb_birth)
                    [~,log_sum_w] = normalizeLogWeights(mbm_upd.w(temp));
                    ppp(end+1,1).w = log_sum_w+log(mbm_upd.track{i}(idx(j)).r);
                    ppp(end,1).xr = mbm_upd.track{i}(idx(j)).xr;
                    ppp(end,1).Cr = mbm_upd.track{i}(idx(j)).Cr;
                    ppp(end,1).V = mbm_upd.track{i}(idx(j)).V;
                    ppp(end,1).v = mbm_upd.track{i}(idx(j)).v;
                    ppp(end,1).alpha = mbm_upd.track{i}(idx(j)).alpha;
                    ppp(end,1).beta = mbm_upd.track{i}(idx(j)).beta;
                    ppp(end,1).s = mbm_upd.track{i}(idx(j)).s;
                    ppp(end,1).t = mbm_upd.track{i}(idx(j)).t;
                end
            end
            %delete these single target hypotheses
            %re-index
            idx_0 = mbm_upd.table(:,i) > 0;
            [idx,~,temp] = unique(mbm_upd.table(idx_0,i));
            mbm_upd.table(idx_0,i) = temp;
            mbm_upd.track{i} = mbm_upd.track{i}(idx);
        end
        %remove empty track
        idx = ~cellfun('isempty',mbm_upd.track);
        mbm_upd.track = mbm_upd.track(idx);
        mbm_upd.table = mbm_upd.table(:,idx);
        if isempty(mbm_upd.table)
            mbm_upd.table = zeros(1,0);
            mbm_upd.track = cell(0,1);
            mbm_upd.w = 0;
        end
        
        %merge rows of global hypothesis look-up table that are same
        if length(mbm_upd.w) > 1
            [mbm_upd.table,~,IC] = unique(mbm_upd.table,'rows');
            n_a = size(mbm_upd.table,1);
            temp = zeros(n_a,1);
            for i = 1:n_a
                [~,temp(i)] = normalizeLogWeights(mbm_upd.w(IC==i));
            end
            mbm_upd.w = temp;
        end
    end
    %% 
    if paras.clutter_est
        nc_hat = gatingBern(mbm_upd,Z{t},measmodel,paras);
        paras.clutter.alpha = paras.clutter.alpha + nc_hat;
        paras.clutter.beta = paras.clutter.beta + 1;
        paras.clutter.lambdaC_hat = paras.clutter.alpha/paras.clutter.beta;
        paras.clutter.c_intensity_hat = measmodel.c_pdf * paras.clutter.lambdaC_hat;
        lambdaC = [lambdaC;paras.clutter.lambdaC_hat];
        Nc_hat = [Nc_hat;nc_hat];
    end
    %%
    dxr = motionmodel.dxr;
    if ~paras.mb_birth
        %remove ppp components with small weights
        idx = ([ppp.w] > paras.pruning.ppp) & ([ppp.alpha]./[ppp.beta] > 1) & ([ppp.v] > 7);
        ppp = ppp(idx);
        if ~paras.mb_approx
            %merge similar ppp components
            ppp = BGGmixtureReduction(ppp,paras.merging.ppp);
        end
    end
    
    %multi-target state estimation
    if paras.estimator == 1
        %find the global hypothesis with the highest weight
        [~,I] = max(mbm_upd.w);
        est{t}.xr = zeros(dxr,0);
        est{t}.X = zeros(2,2,0);
        est{t}.pd = zeros(4,0);
        est{t}.pd_total = (exp([ppp.w])/sum(exp([ppp.w])))*([ppp.s]./([ppp.s]+[ppp.t]))';
        est{t}.label = zeros(1,0);
        est{t}.x_lambda = zeros(1,0);
        if ~isempty(mbm_upd.table)
            hypo_best = mbm_upd.table(I,:);
            for i = 1:length(hypo_best)
                if mbm_upd.table(I,i) > 0
                    %extract state estimate from Bernoulli component with large
                    %enough probability of existence
                    if mbm_upd.track{i}(mbm_upd.table(I,i)).r > paras.estimate.r
                        card_est(t) = card_est(t) + 1;
                        est{t}.xr = [est{t}.xr mbm_upd.track{i}(mbm_upd.table(I,i)).xr];
                        est{t}.X = cat(3,est{t}.X,...
                            mbm_upd.track{i}(mbm_upd.table(I,i)).V...
                            /(mbm_upd.track{i}(mbm_upd.table(I,i)).v-6));
                        pd_hat = mbm_upd.track{i}(mbm_upd.table(I,i)).s/(mbm_upd.track{i}(mbm_upd.table(I,i)).s+mbm_upd.track{i}(mbm_upd.table(I,i)).t);
                        est{t}.pd = [est{t}.pd [mbm_upd.track{i}(mbm_upd.table(I,i)).s;mbm_upd.track{i}(mbm_upd.table(I,i)).t;pd_hat;mbm_upd.track{i}(mbm_upd.table(I,i)).isDetected]];
                        est{t}.label = [est{t}.label mbm_upd.track{i}(mbm_upd.table(I,i)).label];
                        est{t}.x_lambda =[est{t}.x_lambda mbm_upd.track{i}(mbm_upd.table(I,i)).alpha/mbm_upd.track{i}(mbm_upd.table(I,i)).beta];
                    end
                end
            end
            pd_var = est{t}.pd(1,:).*est{t}.pd(2,:)./((est{t}.pd(1,:)+est{t}.pd(2,:)).^2)./(est{t}.pd(1,:)+est{t}.pd(2,:)+1);
            w_var = (1./pd_var)./sum(1./pd_var);
            est{t}.pd_total = w_var*est{t}.pd(3,:)';
            if paras.Pd_est
                measmodel.Pd = est{t}.pd_total;
            end
        end
        est{t}.card = size(est{t}.xr,2);
    elseif paras.estimator == 2
        %compute the cardinality distribution
        [n_a,n_track] = size(mbm_upd.table);
        card_dist = zeros(1,n_track+1);
        pcard = zeros(n_a,n_track+1);
        for i = 1:n_a
            r = [];
            for j = 1:n_track
                if mbm_upd.table(i,j) > 0
                    r = [r mbm_upd.track{j}(mbm_upd.table(i,j)).r];
                end
            end
            if ~isempty(r)
                %avoid numerical underflow
                r(r>1-1e-6) = 1-1e-6;
                pcard(i,1:length(r)+1) = prod(1-r)*poly(-r./(1-r));
            end
             
            card_dist = card_dist + pcard(i,:)*exp(mbm_upd.w(i));
        end
        %obtain the maximum cardinality
        [~,card_max] = max(card_dist);
        %find the global hypothesis with the highest weight and the same
        %MAP cardinality estimate
        [~,a_best] = max(pcard(:,card_max));
        r = zeros(n_track,1);
        xr = zeros(dxr,n_track);
        X = zeros(2,2,n_track);
        for i = 1:n_track
            if mbm_upd.table(a_best,i) > 0
                r(i) = mbm_upd.track{i}(mbm_upd.table(a_best,i)).r;
                xr(:,i) = mbm_upd.track{i}(mbm_upd.table(a_best,i)).xr;
                X(:,:,i) = mbm_upd.track{i}(mbm_upd.table(a_best,i)).V...
                    /(mbm_upd.track{i}(mbm_upd.table(a_best,i)).v-6);
            end
        end
        [~,I] = sort(r,'descend');
        card_est(t) = card_max-1;
        est{t}.xr = xr(:,I(1:card_max-1));
        est{t}.X = X(:,:,I(1:card_max-1));
        est{t}.card = sum(card_dist.*(0:length(card_dist)-1));
    end
    
    %% plot
    % for i = 1:n_track
    %     x = mbm_upd.track{i}(mbm_upd.table(a_best,i)).xr;
    %     handle_extent = plot_extent_iw(x,X,line_style, color, line_width);
    %     H = [H handle_extent];
    % end
    %find the global hypothesis with the highest weight
    if plot_enable_trajectory
        [~,I] = max(mbm_upd.w);
        if ~isempty(mbm_upd.table)
            hypo_best = mbm_upd.table(I,:);
            for i = 1:length(hypo_best)
                if mbm_upd.table(I,i) > 0
                    %extract state estimate from Bernoulli component with large
                    %enough probability of existence
                    if mbm_upd.track{i}(mbm_upd.table(I,i)).r > paras.estimate.r
    
                        x = [est{t}.xr mbm_upd.track{i}(mbm_upd.table(I,i)).xr];
                        X = cat(3,est{t}.X,...
                            mbm_upd.track{i}(mbm_upd.table(I,i)).V...
                            /(mbm_upd.track{i}(mbm_upd.table(I,i)).v-6));
                        for j = 1:size(x,2)
                            handle_extent = plot_extent_iw(x(1:2,j),X(:,:,j),"-", "b", 1);
                            H = [H handle_extent];
                        end
                    end
                end
            end
        end
        pause(0.1);
        delete(H);      
    end
    %%

    %prediction step
    %prediction for detected targets
    mbm.w = mbm_upd.w;
    mbm.table = mbm_upd.table;
    mbm.track = cellfun(@(x) arrayfun(@(x) ...
        BGGbern_pred(x,motionmodel),x),mbm_upd.track,'uniformoutput',false);
    
    if paras.mb_birth
        %add multi-Bernoulli birth
        na = length(mbm.w);
        nt = length(mbm.track);
        mbm.table = [mbm.table ones(na,nb)];
        for i = 1:nb
            mbm.track{nt+i,1}.r = exp(birthmodel(i).w);
            mbm.track{nt+i,1}.xr = birthmodel(i).xr;
            mbm.track{nt+i,1}.Cr = birthmodel(i).Cr;
            mbm.track{nt+i,1}.V = birthmodel(i).V;
            mbm.track{nt+i,1}.v = birthmodel(i).v;
            mbm.track{nt+i,1}.alpha = birthmodel(i).alpha;
            mbm.track{nt+i,1}.beta = birthmodel(i).beta;
            mbm.track{nt+i,1}.s = birthmodel(i).s;
            mbm.track{nt+i,1}.t = birthmodel(i).t;
        end
    else
        %prediction for undetected targets
        ppp = BGGppp_pred(ppp,motionmodel,birthmodel);
    end
    
    t_elapsed(t) = toc;

end

fprintf('\n')

card_est = cellfun(@(x) x.card,est);
card_err = card_est - card;
%evaluate multi-target filtering performance using GOSPAs
d_gospa = zeros(T,1);
decomposed_cost = repmat(struct('localisation',[],'missed',[],'false',[]),T,1);

gt_t = cell(100,1);
for t = 1:T
    x_mat.x = zeros(2,0);
    x_mat.X = zeros(2,2,0);
    for i = 1:length(gt)
        if gt(i).x_bt <= t && gt(i).x_dt >= t
            x_mat.x = [x_mat.x gt(i).xr(1:2,t-gt(i).x_bt+1)];
            x_mat.X = cat(3,x_mat.X,gt(i).X(:,:,t-gt(i).x_bt+1));
        end
    end
    y_mat.x = est{t}.xr;
    y_mat.X = est{t}.X;
    gt_t{t}.x = x_mat.x;
    gt_t{t}.X = x_mat.X;
    [d_gospa(t),~,decomposed_cost(t)] = ...
        GOSPA_extended(x_mat,y_mat,gospa.p,gospa.c,gospa.alpha);
end
% 绘制综合性能曲线
if plot_enable_result
    figure(5);
    clf(5);
    plot(1:T,card,'linewidth',2)
    grid on
    hold on
    plot(1:T,card_est,'linewidth',2)
    xlabel('Time step')
    ylabel('Number of targets')
    legend('True cardinality','Estimated cardinality')
    
    figure(6);
    clf(6);
    subplot(2,2,1)
    plot(1:T,d_gospa,'linewidth',2)
    grid on
    xlabel('Time step')
    ylabel('GOSPA')
    subplot(2,2,2)
    plot(1:T,[decomposed_cost.localisation],'linewidth',2)
    grid on
    xlabel('Time step')
    ylabel('Localisation error')
    subplot(2,2,3)
    plot(1:T,[decomposed_cost.missed],'linewidth',2)
    grid on
    xlabel('Time step')
    ylabel('Misdetection error')
    subplot(2,2,4)
    plot(1:T,[decomposed_cost.false],'linewidth',2)
    grid on
    xlabel('Time step')
    ylabel('False detection error')
end
% 绘制杂波率估计曲线
if paras.clutter_est && plot_enable_result
    figure(7);
    clf(7);
    hold on;grid on;
    plot(1:T,ones(T,1)*measmodel.c_lambda,"k",'DisplayName','Truth');
    plot(1:T,lambdaC,"r","LineWidth",2,'DisplayName','Estimate');
    plot(1:T,Nc_hat,"b");
    xlabel('Time step');
    ylabel('Clutter Rate Estimate');
end
% 绘制检测概率估计曲线
Pd_hat = zeros(T,1);
for t = 1:T
    Pd_hat(t) = est{t}.pd_total;
end
if paras.Pd_est && plot_enable_result
    figure(8);
    clf(8);
    hold on;grid on;
    plot(1:T,ones(T,1)*measmodel.Pd,"k",'DisplayName','Truth');
    plot(1:T,Pd_hat,"r","LineWidth",2,'DisplayName','Estimate');
    xlabel('Time step');
    ylabel('Detection Probability Estimate');
end
%打印结果
fprintf('Mean GOSPA: %.2f\n',mean(d_gospa))
fprintf('Total runtime: %.2f\n',sum(t_elapsed))


est_gt =  repmat(struct('label',0,'x_sur',[],'x_lambda',[],'pd_hat',0,'beta',[],'xr',[],'X',[]),5000,1);
nj = 0;
for t=1:T
    for i = 1:length(est{t}.label)
        idx = est{t}.label(i);
        est_gt(idx).label = idx;
        est_gt(idx).x_sur = [est_gt(idx).x_sur t];
        est_gt(idx).x_lambda = [est_gt(idx).x_lambda est{t}.x_lambda(i)];
        est_gt(idx).pd_hat =  est{t}.pd_total;
        est_gt(idx).beta = [est_gt(idx).beta est{t}.pd(:,i)];
        est_gt(idx).xr = [est_gt(idx).xr est{t}.xr(:,i)];
        est_gt(idx).X = cat(3,est_gt(idx).X, est{t}.X(:,:,i));
    end
end

idx_keep = find([est_gt.label]~=0);
est_gt = est_gt(idx_keep);