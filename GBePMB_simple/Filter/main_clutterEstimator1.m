clc;clear;
rng('default');

addpath('..\clutterEstimate\');
addpath('..\scenarioCreate\');
addpath("..\MEOT-GGIW\");
addpath('..\Third-Party-Code\');
addpath('E:\matlab files\研究生\Coexisting-point-extended-target-PMBM-filter-main\Coexisting-point-extended-target-PMBM-filter-main\Third-party code');


%% motionmodel
%survival probability
Ps =0.99;

%Target kinematic state [x-position,y-position,x-velocity,y-velocity]
%Target extent state [orientation,semi-axis length 1,semi-axis length 2]

%Parameters of a nearly constant velocity motion model
%kinematic state dimension
dxr = 4;
%time interval
Ts = 1;
%transition matrix for kinematic state
Ar = [1 0 Ts 0;
    0 1 0 Ts;
    0 0 1 0;
    0 0 0 1];
%process noise
q = 0.09;
%process noise covariance matrix for kinematic state
Cwr = q*[Ts^3/3 0      Ts^2/2 0;
    0      Ts^3/3 0      Ts^2/2;
    Ts^2/2 0      Ts     0;
    0      Ts^2/2 0      Ts];
%measurement rate parameter used for prediction of gamma distribution
eta = 1.2;
%forgetting factor used for prediction of inverse-Wishart distribution
tau = 20;
% 存活的场景
range_x = [-200 200];
range_y = [-200 200];

%struct representation
motionmodel.Ps = Ps;
motionmodel.Ts = Ts;
motionmodel.dxr = dxr;
motionmodel.Ar = Ar;
motionmodel.Cwr = Cwr;
motionmodel.eta = eta;
motionmodel.tau = tau;
motionmodel.range_x = range_x;
motionmodel.range_y = range_y;

%% birthmodel
 %Poisson birth model with Gaussian mixture intensity
%number of Gaussian components
b_nG = 3;
%Poisson intensit4
birthmodel = repmat(struct('w',log(0.01),'xr',[],'Cr',diag([10 10 4 4])^2,...
    'V',[],'v',100,'alpha',500,'beta',100,'isDetected',false),b_nG,1);
%specify kinematic state (x-position,y-position,x-velocity,y_velocity)
birthmodel(1).xr = [-100 100 0 0]';
birthmodel(2).xr = [-100 -100 0 0]';
birthmodel(3).xr = [100 0 0 0]';

%specify extent state (orientation,two axis lengths)
birthmodel(1).V = (200+6)*diag([1 1])^2;
birthmodel(2).V = (200+6)*diag([1 1])^2;
birthmodel(3).V = (200+6)*diag([1 1])^2;

%Poisson birth rate
b_lambda = sum(exp([birthmodel.w]));

%% measmodel
%detection probability
Pd = 0.9;
%Poisson clutter rate
c_lambda = 60;
%surveillance area
range_x = [-200 200];
range_y = [-200 200];
%uniform distributed clutter density
c_pdf = 1/(range_x(2)-range_x(1))/(range_y(2)-range_y(1));
%Poisson clutter intensity
c_intensity = c_lambda*c_pdf;

%Parameters of a linear Gaussian measurement model
%measurement dimension
dz = 2;
%observation matrix
H = [1 0 0 0;
    0 1 0 0];
%covariance of the multiplicative noise
Ch = diag([1/2 1/2]);
%covariance of the measurement noise
Cv = diag([1/2 1/2]);

%struct representation
measmodel.dz = dz;
measmodel.H = H;
measmodel.Ch = Ch;
measmodel.Cv= Cv;
measmodel.Pd = Pd;
measmodel.c_intensity = c_intensity;
measmodel.c_lambda = c_lambda;
measmodel.c_pdf = c_pdf;


%% groundTruth
T = 100;
% 存活的场景
range_x = [-200 200];
range_y = [-200 200];
%time interval
Ts = 1;
%transition matrix for kinematic state
Ar = [1 0 Ts 0;
0 1 0 Ts;
0 0 1 0;
0 0 0 1];
%create memory to store ground truth
%birth time, death time, target state
gt = struct('x_bt',[],'x_dt',[],'x_lambda',[],'xr',[],'X',[]);

nx = 0;

nx = nx + 1;
gt(nx).x_bt = 1;
gt(nx).x_dt = T;
x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
gt(nx).x_lambda = 5*ones(1,x_dbt);
gt(nx).xr = [-85;85;2.5;-2.5];
gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
for i= 1:x_dbt-1
    x_k1 = Ar*gt(nx).xr(:,end);
    gt(nx).xr = [gt(nx).xr x_k1];
end

nx = nx + 1;
gt(nx).x_bt = 1;
gt(nx).x_dt = 85;
x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
gt(nx).x_lambda = 5*ones(1,x_dbt);
gt(nx).xr = [-80;-80;2.5;2.5];
gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
for i= 1:x_dbt-1
    x_k1 = Ar*gt(nx).xr(:,end);
    gt(nx).xr = [gt(nx).xr x_k1];
end

nx = nx + 1;
gt(nx).x_bt = 50;
gt(nx).x_dt = 70;
x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
gt(nx).x_lambda = 5*ones(1,x_dbt);
gt(nx).xr = [80;20;0;4.2];
gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
for i= 1:x_dbt-1
    x_k1 = Ar*gt(nx).xr(:,end);
    gt(nx).xr = [gt(nx).xr x_k1];
end

nx = nx + 1;
gt(nx).x_bt = 20;
gt(nx).x_dt = 100;
x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
gt(nx).x_lambda = 5*ones(1,x_dbt);
gt(nx).xr = [-120;-120;3;1];
gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
for i= 1:x_dbt-1
    x_k1 = Ar*gt(nx).xr(:,end);
    gt(nx).xr = [gt(nx).xr x_k1];
end

gt_t = repmat(struct('x_lambda',[],'xr',[],'X',[]),T,1);
nt = size(gt,2);
for t = 1:T
    for i = 1:nt
        if gt(i).x_bt<=t && gt(i).x_dt>=t
            gt_t(t).x_lambda = [gt_t(t).x_lambda gt(i).x_lambda(t-gt(i).x_bt+1)];
            gt_t(t).xr = [gt_t(t).xr gt(i).xr(:,t-gt(i).x_bt+1)];
            gt_t(t).X =  cat(3,gt_t(t).X,gt(i).X(:,:,t-gt(i).x_bt+1));
        end
    end
end


%cardinality of multi-target states
card = zeros(T,1);
for i = 1:T
    for j = 1:nt
        if gt(j).x_bt <= i && gt(j).x_dt >= i
            card(i) = card(i) + 1;
        end
    end
end

groundTruth.gt = gt;
groundTruth.gt_time = gt_t; 
groundTruth.card = card;


%% 绘制真实轨迹，绘制量测
plot_enable = 0;
if plot_enable
    plotTrajectory(gt,motionmodel,1);
    % 绘制动态图
    pt = 0.1;%动态图延迟
    plotDynamicTrajectory(gt_t,gt,motionmodel,2,pt);
    % plotDynamicMeas(Z,3,pt);
end

%% 滤波器参数初始化
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
paras.clutter.beta = 10; 
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

%% 预定义空间
N_MC = 100;
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

N_filter = 6;% 滤波器数量
N_MC = 100; % 模拟次数

GOSPA_mean = zeros(N_MC,N_filter);

D_result = repmat(struct('GOSPA',zeros(T,N_MC),...
    'location',zeros(T,N_MC),...
    'miss',zeros(T,N_MC),...
    'false',zeros(T,N_MC),...
    'lambda_c',zeros(T,N_MC),...
    'pd',zeros(T,N_MC),...
    'Card',zeros(T,N_MC), ...
    'Runtime',zeros(1,N_MC))...
    ,N_filter,1);
pd_single = repmat({0.7*ones(T,N_MC)},length(gt),1);

GZ = repmat(struct('gt',[],'gt_time',[],'card',[],'Z',[]),N_MC,1);

%% 产生仿真所需的所有量测
for i = 1:N_MC
    Z = generateMeas(gt,measmodel,motionmodel,T);
    
    GZ(i).gt = gt;
    GZ(i).gt_time = gt_t;
    GZ(i).card = card;
    GZ(i).Z = Z;
    Card_truth(:,i) = card;
end


%%
for i = 1:N_MC
    fprintf('\n%d\n',i);
    
    n_birth = length([ppp.w]);%出生PPP模型数量
    Z = GZ(i).Z;
    
    %% 对比参数设置区 
    array_clutter_est = [1 1 1 0 0 0];
    array_clutter_rate = [20 60 180 20 60 180 ];
    array_clutter_gammascalar = [5 5 5 5 5 5 ];

    array_pd_est = [0 0 0 0 0 0];
    array_pd_rate = [0.9 0.9 0.9 0.9 0.9 0.9];
    array_pd_betascalar = [10 10 10 10 10 10];

    for i_filter = 1:N_filter
        scenarioParas0.clutter_est = array_clutter_est(i_filter);%杂波率估计器开关
        scenarioParas0.alpha_clutter = array_clutter_rate(i_filter)*array_clutter_gammascalar(i_filter);%杂波率估计器gamma分布的初始alpha
        scenarioParas0.beta_clutter = array_clutter_gammascalar(i_filter);%杂波率估计器gamma分布的初始beta
        % 杂波率初始值
        scenarioParas0.lambdaC_hat = scenarioParas0.alpha_clutter/scenarioParas0.beta_clutter;
        scenarioParas0.Pd_est = array_pd_est(i_filter);%检测概率估计器开关
        scenarioParas0.s_Pd =array_pd_betascalar(i_filter)*array_pd_rate(i_filter)*ones(n_birth,1);%检测概率估计器beta分布的初始s
        scenarioParas0.t_Pd = array_pd_betascalar(i_filter)*(1-array_pd_rate(i_filter))*ones(n_birth,1);%检测概率估计器beta分布的初始t 
        w = exp([birthmodel.w])/sum(exp([birthmodel.w]));
        %检测概率初始值pl
        scenarioParas0.Pd_hat = w*(scenarioParas0.s_Pd./(scenarioParas0.t_Pd+scenarioParas0.s_Pd));
        % 场景参数初始化
        [paras,birthmodel,ppp] = scenarioParasInit(scenarioParas0,paras,measmodel,birthmodel,ppp);
        % 滤波算法
        [d_gospa,decomposed_cost,card_est,lambdaC,pd_hat,est,est_gt,runtime] = ...
        filter2(birthmodel,measmodel,motionmodel,ppp,mbm,paras,groundTruth,Z);
        % 储存结果
        D_result(i_filter).GOSPA(:,i) = d_gospa;
        D_result(i_filter).location(:,i) = [decomposed_cost.localisation]';
        D_result(i_filter).miss(:,i) = [decomposed_cost.missed]';
        D_result(i_filter).false(:,i) = [decomposed_cost.false]';
        D_result(i_filter).lambda_c(:,i) = lambdaC;
        D_result(i_filter).pd(:,i) = pd_hat;
        D_result(i_filter).Card(:,i) = card_est;
        D_result(i_filter).runtime(:,i) = runtime;
        GOSPA_mean(i,i_filter) = mean(d_gospa,1);
    end

    %%  绘图区
    % /===============坐标标签=====================
    x_tick = [];
    if any(array_clutter_est) && ~any(array_pd_est)
        for i_filter = 1:N_filter
            if array_clutter_est(i_filter)
                x_tick = [x_tick strcat("Est-init:",string(array_clutter_rate(i_filter)))];
            else
                x_tick = [x_tick strcat("init:",string(array_clutter_rate(i_filter)))];
            end
        end
    elseif any(array_pd_est) && ~any(array_clutter_est)
        for i_filter = 1:N_filter
            if array_pd_est(i_filter)
                x_tick = [x_tick strcat("Est-init:",string(array_pd_rate(i_filter)))];
            else
                x_tick = [x_tick strcat("init:",string(array_pd_rate(i_filter)))];
            end
        end
    elseif any(array_pd_est) && any(array_clutter_est)
        for i_filter = 1:N_filter
            if array_pd_est(i_filter)
                x_tick = [x_tick strcat("Est-init:",string(array_clutter_rate(i_filter)),"/",string(array_pd_rate(i_filter)))];
            else
                x_tick = [x_tick strcat("init:",string(array_clutter_rate(i_filter)),"/",string(array_pd_rate(i_filter)))];
            end
        end
    else
        for i_filter = 1:N_filter
            x_tick = [x_tick strcat("init:",string(array_clutter_rate(i_filter)),"/",string(array_pd_rate(i_filter)))];
        end
    end
    %===============坐标标签=====================/

    % /=============GOSPA==============
    fig_id = 10;
    fig_id  = fig_id+1;figure(fig_id);clf(fig_id);hold on;box on;grid on;legend;
    tpye_set = ["r-*","g-*","b-*","c-^","k-^","m-^"];
    for i_filter = 1:N_filter
        plot(1:T,sum(D_result(i_filter).GOSPA,2)./i, ...
            tpye_set(i_filter),...
            "LineWidth",1, ...
            'DisplayName',x_tick(i_filter));
    end
    xlabel('Time step');
    ylabel('GOSPA');
    % =============GOSPA=================/
    
    % /===============绘制目标基数估计图================
    fig_id  = fig_id+1;figure(fig_id);clf(fig_id);hold on;box on;grid on;legend;
    plot(1:T,mean(Card_truth,2),"k--","LineWidth",2,'DisplayName','Truth');
    
    for i_filter = 1:N_filter
        plot(1:T,sum(D_result(i_filter).Card,2)./i, ...
            tpye_set(i_filter),...
            "LineWidth",1, ...
            'DisplayName',x_tick(i_filter));
    end
    xlabel('Time step');
    ylabel('目标基数');
    % ===============绘制目标基数估计图================/
    
    % /==============绘制杂波率估计图=================
    fig_id  = fig_id+1;figure(fig_id);clf(fig_id);hold on;box on;grid on;legend;
    plot(1:T,ones(T,1)*measmodel.c_lambda,"k","LineWidth",1,'DisplayName','Truth');
    plot(1:T,sum(D_result(1).lambda_c,2)./i,"r","LineWidth",2,'DisplayName',x_tick(1));
    plot(1:T,sum(D_result(2).lambda_c,2)./i,"g","LineWidth",2,'DisplayName',x_tick(2));
    plot(1:T,sum(D_result(3).lambda_c,2)./i,"b","LineWidth",2,'DisplayName',x_tick(3));
    xlabel('Time step');
    ylabel('Clutter Rate Estimate');
    % ==============绘制杂波率估计图=================/

    %/==============绘制检测概率估计图===============
    fig_id  = fig_id+1;figure(fig_id);clf(fig_id);hold on;box on;grid on;legend;
    plot(1:T,ones(T,1)*measmodel.Pd,"k","LineWidth",2,'DisplayName','Truth');
    plot(1:T,sum(D_result(1).pd,2)./i,"r","LineWidth",2,'DisplayName',x_tick(1));
    plot(1:T,sum(D_result(2).pd,2)./i,"g","LineWidth",2,'DisplayName',x_tick(2));
    plot(1:T,sum(D_result(3).pd,2)./i,"b","LineWidth",2,'DisplayName',x_tick(3));
    xlabel('Time step');
    ylabel('Detection Probability Estimate');
    % ==============绘制检测概率估计图===============/
end

%% 对比结果

aa=[mean(mean(D_result(1).GOSPA,2),1);
mean(mean(D_result(2).GOSPA,2),1);
mean(mean(D_result(3).GOSPA,2),1);
mean(mean(D_result(4).GOSPA,2),1);
mean(mean(D_result(5).GOSPA,2),1);
mean(mean(D_result(6).GOSPA,2),1);
];
fig_id  = fig_id+1;figure(fig_id);clf(fig_id);hold on;box on;grid on;

x = categorical(x_tick);
x = reordercats(x,x_tick);
hb  =bar(x,aa);
hb.FaceColor = 'flat';
hb.CData(1,:)=[1,0,0]';
hb.CData(2,:)=[1,0,0]';
hb.CData(4,:)=[1,0.5,0]';
ylim([min(aa)-0.02*min(aa),max(aa)+0.02*min(aa)]);
ylabel("平均GOSP0A");

fprintf('Mean GOSPA improve:\n %.2f\n%.2f\n%.2f\n%.2f\n%.2f\n',(aa-aa(4))*100/aa(4));

fig_id  = fig_id+1;figure(fig_id);clf(fig_id);hold on;box on;grid on;
x = categorical( [x_tick,"Imp-"+x_tick(1:2)]);
x = reordercats(x,[x_tick,"Imp-"+x_tick(1:2)]);
hb  =bar(x,[(aa-aa(4))*100/aa(4);[(aa(3)-aa(1))*100/aa(1);(aa(5)-aa(2))*100/aa(2)]]);
hb.FaceColor = 'flat';
hb.CData(1,:)=[1,0,0]';
hb.CData(2,:)=[1,0,0]';
hb.CData(4,:)=[1,0.5,0]';
hb.CData(6,:)=[0,1,0]';
hb.CData(7,:)=[0,1,0]';
ylabel("平均GOSP0A improve");


fprintf('compare GOSPA improve:\n %.2f\n%.2f\n',[(aa(3)-aa(1))*100/aa(1),(aa(5)-aa(2))*100/aa(2)]);
fprintf("\nGOSPA:\n")
fprintf('%f\n',aa);