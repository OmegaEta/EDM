clc;clear;dbstop if error;
% clc;clear;close all;dbstop if error;
dbstop if error ;
addpath(genpath('.\Scenario'));addpath(genpath('.\lib'));
addpath(genpath('.\lib_EDM_PMBM'));addpath(genpath('.\lib_EDM'))
scenario = 2;
%Parameter setting.
if scenario == 1
    modelparas1;
elseif scenario == 2
    modelparas2;
end

%Number of Monte Carlo Simulations
% numMC = length(Scenario.Z);
numMC=100;
fig = 1;
%Parameters used in GOSPA metric
c = 20;
p = 1;

%Number of time steps
K = model.K;

GOSPA = zeros(K,5,numMC);
GOSPA_ = zeros(K,5,numMC);
trajectoryEstimates = cell(numMC,1);
simulation_time = zeros(numMC,2);

%% EDM-PMBM
if 1
model.birth.GGIW = model.birth.irregularGGIW;
for t = 1:numMC
    Z = Scenario.Z{t};
    % Initialisation
    PPP.w = log(model.birth.w);
    PPP.GGIW = model.birth.GGIW;
    MBM.w = [];     % Global hypotheses weights
    MBM.track = {}; % Locl hypotheses trees
    MBM.table = []; % Global hypotheses look-up table
     figure(1);
        clf(1);
    estimates = cell(K,1);
    
    tic
    for k = 1:K
        %Print info
        pause(0);

        axis(model.Scenario_range);
        
        hold on;
        disp(['EDM-PMBM' 'Scenario' + string(scenario)  'Monte Carlo runs' + string(t) 'TimeStep' + string(k)]);
        
        %Update step
        [PPP,MBM] = updatePMBM(PPP,MBM,Z{k},k,model);
        
        %Extract estimates (both estimate of the current time and the
        %estimate of the full trajectory) 
        [estimates_irregular{k},trajectoryEstimates_irregular{t}{k}] = estimator_irregular(MBM,model);
        
%         save estimates_irregular_forTruth3.mat estimates_irregular
        %Evaluate filtering performance using GOSPA
        GOSPA(k,:,t) = GOSPAmetric_EDM(estimates_irregular{k},groundTruth_target_irregular{k},c,p);
if fig == 1 && t==1
        figure(1);
         if ~isempty(Z{k})
            plot(Z{k}(1,:),Z{k}(2,:),"k.",'HandleVisibility','off');
         end
%         for I=1:length(MBM.track)
%             if(~isempty(Z{k}))
%                 for J=1:length(MBM.track{I})
%                     plot(Z{k}(1,MBM.track{I}(J).assocHistory(end).meas(:)),Z{k}(2,MBM.track{I}(J).assocHistory(end).meas(:)),"k.",'HandleVisibility','off');
% %                     for k1 = 1:length(MBM.track{I})
% %                         plot(MBM.track{I}(k1).Bern.GGIW(end).m(1),MBM.track{I}(k1).Bern.GGIW(end).m(2),'.');
% %                     end
%                 end
%             end
%         end
        i0=size(estimates_irregular{k}.g,2);
        for j = 1:i0
            figure(1)
            hold on
%              [x, y] = Sigmacircle(estimates{k}.x(1,j),estimates{k}.x(2,j),estimates{k}.X(:,:,j),j+4-1,j+3);
             [x_ir, y_ir] = irregularcircle({estimates_irregular{k}.r{j,1}},estimates_irregular{k}.x(:,j),model.direction_angle,j);

%              [x_i, y_i] = irregularcircle({groundTruth_target_irregular{k}.r{j,1}},groundTruth_target_irregular{k}.x(:,j),model.direction_angle,7);
        end
end
         drawnow;

        if k < K
            [PPP,MBM] = predictPMBM(PPP,MBM,model);
%             ys = evaluateFourierSeries(MBM.track{2, 1}.Bern.GGIW(end).Shape_coefficients.a,MBM.track{2, 1}.Bern.GGIW(end).Shape_coefficients.b,model.direction_angle);
%             y = evaluateFourierSeries(MBM.track{1, 1}.Bern.GGIW(end).shape.dilation_coefficients.a,MBM.track{1, 1}.Bern.GGIW(end).shape.dilation_coefficients.b,model.direction_angle);
        end
    end
    simulation_time(t,1) = toc;
end
GOSPA02_EDMPMBM = sum(GOSPA,3) / numMC;
for i=1:4
    figure(1+i)
    plot(1:K,GOSPA02_EDM(:,i),'-r',LineWidth=2);
    hold on;
    title(model.GOSPAtitles{i});%
end
rmpath(genpath('.\lib_EDM_PMBM'));
end

%% GGIW-PMBM
if 1
addpath(genpath('.\lib_GGIW_PMBM'));addpath(genpath('.\lib'));
for t = 1:numMC
    Z = Scenario.Z{t};
    % Initialisation
    PPP.w = log(model.birth.w);
    PPP.GGIW = model.birth.GGIW;
    MBM.w = [];     % Global hypotheses weights
    MBM.track = {}; % Locl hypotheses trees
    MBM.table = []; % Global hypotheses look-up table
    
    estimates = cell(K,1);
    
    tic
    for k = 1:K
        %Print info
        pause(0);
        figure(1);
%         clf(1);
        axis(model.Scenario_range);
        
        hold on;
         disp(['GGIW-PMBM' 'Scenario' + string(scenario)  'Monte Carlo runs' + string(t) 'TimeStep' + string(k)]);
 
        %Update step
        [PPP,MBM] = updatePMBM(PPP,MBM,Z{k},k,model);
                
        %Extract estimates (both estimate of the current time and the
        %estimate of the full trajectory) 
        [estimates{k},trajectoryEstimates{t}{k}] = estimator(MBM,model);

        %Evaluate filtering performance using GOSPA
        GOSPA_(k,:,t) = GOSPAmetric_EDM(estimates{k},groundTruth_target_irregular{k},c,p);

if fig == 1 && t==1
        figure(1);
         if ~isempty(Z{k})
            plot(Z{k}(1,:),Z{k}(2,:),"k.",'HandleVisibility','off');
         end

        i0=size(estimates{k}.g,2);
        for j = 1:i0
            figure(1);
            hold on
            [x, y] = Sigmacircle(estimates{k}.x(1,j),estimates{k}.x(2,j),estimates{k}.X(:,:,j),2,j);
%             [x,y] = irregularcircle({estimates{k}.r{j,1}},estimates{k}.x(:,j),model.direction_angle,j);

        end
end
         drawnow;
        %Prediction Step 
        if k < K
            [PPP,MBM] = predictPMBM(PPP,MBM,model);
        end
    end
    simulation_time(t,2) = toc;
end
GOSPA02_GGIW = sum(GOSPA_,3) / numMC;

for i=1:4
    figure(1+i)
    plot(1:K,GOSPA02_GGIWPMBM(:,i),'-b');
    hold on;
    title(model.GOSPAtitles{i}); 
    legend('GGIW-PMBM'); 
end

rmpath(genpath('.\lib_PMBM'));
end

%% KS-PMBM
if 1
addpath(genpath('.\ks_pmbm'));

%Parameter setting.
if scenario == 1
    modelparas1;
elseif scenario == 2
    modelparas2;
end

GOSPA_improved = zeros(K,5,numMC);
% GOSPA_original = zeros(K,4,numMC);
trajectoryEstimates_improved = cell(numMC,1);
% trajectoryEstimates_original = cell(numMC,1);
simulation_time = zeros(numMC,1);

for t = 1:numMC
    % Initialisation
    Z = Scenario.Z{t};

    PPP_improved.w = log(model.birth.w);
    PPP_improved.GGIW = model.birth.GGIW;

    MBM_improved.w = [];     % Global hypotheses weights
    MBM_improved.track = {}; % Locl hypotheses trees
    MBM_improved.table = []; % Global hypotheses look-up table
 
    estimates_improved = cell(K,1);
    estimates = cell(K,1);
    trajectoryEstimates_improved{t} = cell(K,1);

    cost_time = zeros(K,1);
    tic
    for k = 1:K
        disp(['KS-PMBM' 'Scenario' + string(scenario)  'Monte Carlo runs' + string(t) 'TimeStep' + string(k)]);
        tic;
        pause(0);
        figure(1);
        axis(model.Scenario_range);
        
        %Update step
        [PPP_improved,MBM_improved] = updatePMBM_improved(PPP_improved,MBM_improved,Z{k},k,model,scenario);%改进算法

        %Extract estimates (both estimate of the current time and the
        %estimate of the full trajectory) 
        [estimates_ks{k},trajectoryEstimates_improved{t}{k}] = estimator(MBM_improved,model);

%         GOSPA_improved(k,:,t) = GOSPAmetric_EDM(estimates_improved{k},groundTruth_target_irregular{k},c,p);
        GOSPA_ks(k,:,t) = GOSPAmetric_EDM(estimates_ks{k},groundTruth_target_irregular{k},c,p);
        
        % 画出跟踪结果
if fig == 1 && t==1
        figure(1);
         if ~isempty(Z{k})
            plot(Z{k}(1,:),Z{k}(2,:),"k.",'HandleVisibility','off');
         end

         estimates_improved = estimates_ks;
        i0=size(estimates_improved{k}.g,2);
        for j = 1:i0
            figure(1);
            hold on
            [x, y] = Sigmacircle(estimates_improved{k}.x(1,j),estimates_improved{k}.x(2,j),estimates_improved{k}.X(:,:,j),2,j);
%             [x,y] = irregularcircle({estimates{k}.r{j,1}},estimates{k}.x(:,j),model.direction_angle,j);
        end
end
         %Prediction Step
        if k < K
            [PPP_improved,MBM_improved] = predictPMBM(PPP_improved,MBM_improved,model);
        end

        cost_time(k) = toc;
    end
    
    simulation_time(t) = toc;
end

 %100个时刻的GOSPA的平均
GOSPA02_ks = sum(GOSPA_ks,3) / numMC;
for i=1:4
 figure(1+i)
 plot(1:K,GOSPA02_ks(:,i),'-g',LineWidth=1);
 hold on;
 title(model.GOSPAtitles{i}); 
end

rmpath(genpath('.\ks_pmbm'));
end


%% BGGIW-PMB and GBePMB
if 1
addpath('.\GBePMB_simple\Filter');
addpath('.\lib');
addpath('.\GBePMB_simple\clutterEstimate\');
addpath('.\GBePMB_simple\scenarioCreate\');
addpath(".\GBePMB_simple\MEOT-GGIW\");
addpath('.\GBePMB_simple\Third-Party-Code\');

% 模型参数初始化
N_MC = numMC;
motionmodel = paraOfMotionmodel(1);
birthmodel = paraOfBirthmodel(7);
measmodel = paraOfMeasmodel(2);

SC=scenario;

if SC==1
T=50;
modelparas1_GB;

motionmodel = paraOfMotionmodel(1);
birthmodel = paraOfBirthmodel(6);
measmodel = paraOfMeasmodel(1);
else
T=100;
modelparas2_GB;

motionmodel = paraOfMotionmodel(2);
birthmodel = paraOfBirthmodel(7);
measmodel = paraOfMeasmodel(2);
end

% 产生真实轨迹
[gt,gt_time,card] = generateGroundtruth(7,birthmodel,motionmodel,T);
groundTruth.gt = gt;
groundTruth.gt_time = gt_time; 
groundTruth.card = card;
Z = Scenario.Z;

% 滤波器参数初始化
[paras,ppp,mbm] = parasSetting(birthmodel,measmodel,1);

% 预定义空间
N_filter = 3;
D_gospa1 = zeros(T,N_MC);
D_gospa2 = zeros(T,N_MC);
D_gospa3 = zeros(T,N_MC);
D_gospa4 = zeros(T,N_MC);
% location/missed/false error
Decomposed_cost1 = cell(N_MC,1);
Decomposed_cost2 = cell(N_MC,1);
Decomposed_cost3 = cell(N_MC,1);
Decomposed_cost4 = cell(N_MC,1);
% 目标基数
Card_est1 = zeros(T,N_MC);
Card_est2 = zeros(T,N_MC);
Card_est3 = zeros(T,N_MC);
Card_est4 = zeros(T,N_MC);
Card_truth = zeros(T,N_MC);%真实目标基数

Lambda_C = zeros(T,N_MC);%杂波估计结果
Pd_hat = zeros(T,N_MC);%检测概率估计结果

GOSPA_mean = zeros(N_MC,N_filter);

D_result = repmat(struct('GOSPA',zeros(T,N_MC),...
    'location',zeros(T,N_MC),...
    'miss',zeros(T,N_MC),...
    'false',zeros(T,N_MC),...
    'lambda_c',zeros(T,N_MC),...
    'pd',zeros(T,N_MC),...
    'Card',zeros(T,N_MC),...
    'runtime',zeros(1,N_MC), ...
    'est_gt',{cell(N_MC,1)})...
    ,N_filter,1);
pd_single = repmat({0.7*ones(T,N_MC)},length(gt),1);


dynamicGenerate = 1;
if dynamicGenerate==0
    load('GZ60_07.mat');
end

GOSPA_ = zeros(T,5,N_MC);

% 对比滤波器

for i_filter = 1:N_filter
    if i_filter==2
       continue
    end

    n_birth = length([ppp.w]); %出生PPP模型数量
    %%对比参数设置区
    array_clutter_est = [0 0 0];
    array_clutter_rate = [60 60 60];
    array_clutter_gammascalar = [5 5 5];

    array_pd_est = [1 1 0];
    array_pd_rate = [0.6 0.6 0.6];
    array_pd_betascalar = [10 10 10];

    for i = 1:N_MC
    fprintf('\n%d\n',i);
    % plotTrajectory_addition(gt,motionmodel,1);

        scenarioParas0.clutter_est = array_clutter_est(i_filter);%杂波率估计器开关
        scenarioParas0.alpha_clutter = array_clutter_rate(i_filter)*array_clutter_gammascalar(i_filter);%杂波率估计器gamma分布的初始alpha
        scenarioParas0.beta_clutter = array_clutter_gammascalar(i_filter);%杂波率估计器gamma分布的初始beta
        % 杂波率初始值
        scenarioParas0.lambdaC_hat = scenarioParas0.alpha_clutter/scenarioParas0.beta_clutter;
        scenarioParas0.Pd_est = array_pd_est(i_filter);%检测概率估计器开关
        scenarioParas0.s_Pd =array_pd_betascalar(i_filter)*array_pd_rate(i_filter)*ones(n_birth,1);%检测概率估计器beta分布的初始s
        scenarioParas0.t_Pd = array_pd_betascalar(i_filter)*(1-array_pd_rate(i_filter))*ones(n_birth,1);%检测概率估计器beta分布的初始t 
        w = exp([birthmodel.w])/sum(exp([birthmodel.w]));
        %检测概率初始值
        scenarioParas0.Pd_hat = w*(scenarioParas0.s_Pd./(scenarioParas0.t_Pd+scenarioParas0.s_Pd));
        % 场景参数初始化
        [paras,birthmodel,~,ppp] = scenarioParasInit(scenarioParas0,paras,measmodel,birthmodel,ppp);
        if i_filter==1
            paras.Pd_estimator = 1;%第一个滤波器使用BGGIW更新
        else
            paras.Pd_estimator = 2;%其余滤波器使用GBePMB更新
        end

        % 滤波算法
        [d_gospa,decomposed_cost,card_est,lambdaC,pd_hat,est,est_gt,runtime,GOSPA] = ...
        filter2(birthmodel,measmodel,motionmodel,ppp,mbm,paras,groundTruth,Z{N_MC},model,groundTruth_target_irregular,T);
        GOSPA_(:,:,i) = GOSPA(:,:);

        [SC i_filter i]
        if i_filter==1
            disp(['BGGIW-PMB' 'Scenario' + string(SC)  'Monte Carlo runs' + string(i)]);
        else
            disp(['GBePMB' 'Scenario' + string(SC)  'Monte Carlo runs' + string(i)]);
        end
    end

    if i_filter==1
         GOSPA02_BGGIW = sum(GOSPA_,3) / N_MC;
         GOSPA02_ = GOSPA02_BGGIW;
    else
         GOSPA02_GBePMB = sum(GOSPA_,3) / N_MC;
         GOSPA02_ = GOSPA02_GBePMB;
    end

    tpye_set = ["m--","b--","G:"];
    for i=1:4
     figure(1+i)
     plot(1:T,GOSPA02_(:,i),tpye_set(i_filter),LineWidth=2);
     hold on;
     title(model.GOSPAtitles{i});
    end

end

rmpath('.\GBePMB_simple\Filter');
rmpath('.\GBePMB_simple\clutterEstimate\');
rmpath('.\GBePMB_simple\scenarioCreate\');
rmpath(".\GBePMB_simple\MEOT-GGIW\");
rmpath('.\GBePMB_simple\Third-Party-Code\');

% 创建一个包含字符串和数字的元胞数组
end

