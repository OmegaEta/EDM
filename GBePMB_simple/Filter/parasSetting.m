function [paras,ppp,mbm] = parasSetting(birthmodel,measmodel,scenario)
% 滤波的参数设置
if scenario==1
    %gate size in probability
    paras.gating.Pg = 0.999;
    paras.gating.size = chi2inv(paras.gating.Pg,measmodel.dz);
    %hyperparameters in DBSCAN
    paras.dbscan.max_dist = 5;
    paras.dbscan.min_dist = 0.1;
    paras.dbscan.grid_dist = 0.1;
    %data association algorithm
    %method 1: Murty's algorithm
    %method 2: Gibbs sampling
    paras.data_assoc = 1;
    %number of iterations in Gibbs sampling
    paras.gibbs = 100;
    
    %whether to use multi-Bernoulli birth model (the default is ppp birth)
    paras.mb_birth = false;
    paras.ppp_birth = ~paras.mb_birth;
    
    %pruning threshold for global hypotheses
    paras.pruning.w = log(1e-2);
    %pruning threshold for ppp intensity
    paras.pruning.ppp = log(1e-3);
    %merging threshold for ppp intensity
    paras.merging.ppp = 2;
    %cap of global hypotheses
    paras.cap.w = 100;
    %pruning threshold for Bernoulli
    paras.pruning.r = 1e-3;
    %whether to perform recycling
    paras.recycle = true;
    if paras.mb_birth
        paras.recycle = false;
    end
    if paras.recycle
        paras.pruning.r = 1e-1;
    end
    
    %whether to perform MB approximation
    paras.mb_approx = true;
    %M-best assignments in Murty
    paras.M = 20;
    %MB approximation methods
    %method 1: track-oriented merging
    %method 2: SJPDA type merging
    %method 3: LP merging
    paras.mb_merge = 1;
    
    %two different ways to initiate new tracks
    %1: initiate a new track for each measurement
    %2: initiate a new track for each cluster (this one is only valid for
    %track-oriented merging)
    paras.new_track_init = 1;
    
    %convergence threshold for variational approximation
    paras.vb_threshold = 1e-2;
    
    %estimator to extract multi-target state
    %estimator 1: extract state from global hypothesis with the highest weight
    %and Bernoulli components with large enough probability of existence
    %estimator 2: MAP cardinality estimator
    %threshold to extract state estimate from Bernoulli
    paras.estimator = 1;
    paras.estimate.r = 0.5;
    
    %杂波率的gamma估计
    paras.clutter.alpha = 600;
    paras.clutter.beta = 5; 
    paras.clutter.lambdaC_hat = paras.clutter.alpha/paras.clutter.beta; 
    paras.clutter.c_intensity_hat = paras.clutter.lambdaC_hat*measmodel.c_pdf;
    paras.clutter_est = true;
    %检查概率的beta估计
    paras.Pd_est = 0;

    %initialisation parameters
    if paras.ppp_birth
        %global hypothesis weight in logarithm
        mbm.w = 0;
        %global hypothesis look-up table
        mbm.table = zeros(1,0);
        %local hypothesis trees (collections of single target hypotheses)
        mbm.track = cell(0,1);
        %each Bernoulli is parameterised by 1) existence probability r, 2) mean and
        %covariance of the kinematic state xr, Cr, 3) parameters of the
        %extent state V, v, 4) parameters of gamma distribution alpha, beta.
        %PPP for undetected targets, initialised using birth model
        ppp = birthmodel;
    else
        %initial setting for multi-Bernoulli birth model
        mbm.w = 0;
        nb = length(birthmodel);
        mbm.table = ones(1,nb);
        mbm.track = cell(nb,1);
        for i = 1:nb
            mbm.track{i}.r = exp(birthmodel(i).w);
            mbm.track{i}.xr = birthmodel(i).xr;
            mbm.track{i}.Cr = birthmodel(i).Cr;
            mbm.track{i}.V = birthmodel(i).V;
            mbm.track{i}.v = birthmodel(i).v;
            mbm.track{i}.alpha = birthmodel(i).alpha;
            mbm.track{i}.beta = birthmodel(i).beta;
        end
    end

elseif scenario==2
        %gate size in probability
    paras.gating.Pg = 0.999;
    paras.gating.size = chi2inv(paras.gating.Pg,measmodel.dz);
    %hyperparameters in DBSCAN
    paras.dbscan.max_dist = 5;
    paras.dbscan.min_dist = 0.1;
    paras.dbscan.grid_dist = 0.1;
    %data association algorithm
    %method 1: Murty's algorithm
    %method 2: Gibbs sampling
    paras.data_assoc = 1;
    %number of iterations in Gibbs sampling
    paras.gibbs = 100;
    
    %whether to use multi-Bernoulli birth model (the default is ppp birth)
    paras.mb_birth = false;
    paras.ppp_birth = ~paras.mb_birth;
    
    %pruning threshold for global hypotheses
    paras.pruning.w = log(1e-2);
    %pruning threshold for ppp intensity
    paras.pruning.ppp = log(1e-5); %-3
    %merging threshold for ppp intensity
    paras.merging.ppp = 2;
    %cap of global hypotheses
    paras.cap.w = 100;
    %pruning threshold for Bernoulli
    paras.pruning.r = 1e-3;
    %whether to perform recycling
    paras.recycle = true;
    if paras.mb_birth
        paras.recycle = false;
    end
    if paras.recycle
        paras.pruning.r = 1e-1;
    end
    
    %whether to perform MB approximation
    paras.mb_approx = true;
    %M-best assignments in Murty
    paras.M = 20;
    %MB approximation methods
    %method 1: track-oriented merging
    %method 2: SJPDA type merging
    %method 3: LP merging
    paras.mb_merge = 1;
    
    %two different ways to initiate new tracks
    %1: initiate a new track for each measurement
    %2: initiate a new track for each cluster (this one is only valid for
    %track-oriented merging)
    paras.new_track_init = 1;
    
    %convergence threshold for variational approximation
    paras.vb_threshold = 1e-2;
    
    %estimator to extract multi-target state
    %estimator 1: extract state from global hypothesis with the highest weight
    %and Bernoulli components with large enough probability of existence
    %estimator 2: MAP cardinality estimator
    %threshold to extract state estimate from Bernoulli
    paras.estimator = 1;
    paras.estimate.r = 0.5;
    
    %杂波率的gamma估计
    paras.clutter.alpha = 600;
    paras.clutter.beta = 10; 
    paras.clutter.lambdaC_hat = paras.clutter.alpha/paras.clutter.beta; 
    paras.clutter.c_intensity_hat = paras.clutter.lambdaC_hat*measmodel.c_pdf;
    paras.clutter_est = false;

    %检查概率的beta估计
    paras.Pd_est = 1;

    %initialisation parameters
    if paras.ppp_birth
        %global hypothesis weight in logarithm
        mbm.w = 0;
        %global hypothesis look-up table
        mbm.table = zeros(1,0);
        %local hypothesis trees (collections of single target hypotheses)
        mbm.track = cell(0,1);
        %each Bernoulli is parameterised by 1) existence probability r, 2) mean and
        %covariance of the kinematic state xr, Cr, 3) parameters of the
        %extent state V, v, 4) parameters of gamma distribution alpha, beta.
        %PPP for undetected targets, initialised using birth model
        ppp = birthmodel;
    else
        %initial setting for multi-Bernoulli birth model
        mbm.w = 0;
        nb = length(birthmodel);
        mbm.table = ones(1,nb);
        mbm.track = cell(nb,1);
        for i = 1:nb
            mbm.track{i}.r = exp(birthmodel(i).w);
            mbm.track{i}.xr = birthmodel(i).xr;
            mbm.track{i}.Cr = birthmodel(i).Cr;
            mbm.track{i}.V = birthmodel(i).V;
            mbm.track{i}.v = birthmodel(i).v;
            mbm.track{i}.alpha = birthmodel(i).alpha;
            mbm.track{i}.beta = birthmodel(i).beta;
        end
    end
end
end

