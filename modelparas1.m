load('Scenario\Scenario_InterSect_Irregular1.mat');
load('Scenario\Scenario_Intersect_TruthIrregular.mat')
model.Scenario_range = [-180,180,-100,100];
K = length(Scenario.Z{1});
model.K = K;

% Effective window length for the gamma prediction
w_e_gamma = 20;
% Effective window length for the extent prediction
w_e_extent = 15;

model.tao = 1/(log(w_e_extent)-log(w_e_extent-1));
model.eta = 1/(1-1/w_e_gamma);
model.Ts = 1;   %sampling interval
sigma_v = 0.1;  %standard deviation of motion noise
sigma_r = 0.1;  %standard deviation of measurement noise
model.motionmodel = motionmodel.cvmodel(model.Ts,sigma_v);
model.measmodel = measmodel.cvmeasmodel(sigma_r);

%target existence probability
model.Ps = 0.99;
%target detection probability
model.Pd = Scenario.detection_prob;

%range of the surveillance area
range_c = [-1 1;-1 1]*200;
%Poisson false alarm (clutter) rate
% lambda_c = Scenario.false_alarm_rate;
lambda_c = 60;
%Poisson clutter intensity
model.lambda_fa = lambda_c/prod(range_c(:,2)-range_c(:,1));

% Gating parameters
Pg = 0.999;
model.gamma= chi2inv(Pg,model.measmodel.d);
model.Qd = 1 - model.Pd*Pg;

% Thresholds
model.threshold_r = 1e-2;   %existence probability of Bernoulli component
model.threshold_u = 1e-2;   %weight of mixture component in PPP
model.threshold_w = 0.05;   %1e-2 weight of global hypothesis (multi-Bernoulli)
model.threshold_s = 1e-4;   %weight of the trajectory is still alive
model.recycle = 1e-1;       %recycling threshold
model.merge = 4;            %merge threshold used to merge similar GGIWs
model.M = 100;              %cap of number of MBM components in PMBM
model.num_iterations = 3;   %controls the number of iterations used in SO
model.max_repetition = 1;   %controls the number of iterations used in SO

%extract target state from Bernoulli components with existence probability
%larger than this threshold 
model.exist_r = 0.5;        

% target initial state
nbirths = 4;
model.nbirths = nbirths;
xstart = zeros(model.motionmodel.d,nbirths);

xstart(:,1) = [-192 75 0 0];
xstart(:,2) = [-192 10 0 0];
xstart(:,3) = [-192 -30 0 0];
xstart(:,4) = [500 500 0 0];

%Birth model
d = 2;
model.direction = 36;
model.direction_angle = linspace(0, 2*pi, model.direction);
%% EDM参数
model.coefficients_dilation.a = [1 zeros(1,model.direction)];model.coefficients_dilation.a(end) = [];
model.coefficients_dilation.b = zeros(1,model.direction+1);model.coefficients_dilation.b(end) = [];

F_a = 8; F_a = repmat(F_a,1,model.direction);
F_b = 6; F_b = repmat(F_b,1,model.direction);
F_P = diag([4,4,4,4]);%shape_p
F_v_pm = 26.5; F_v = repmat(F_v_pm,1,model.direction);
F_V_pm = 55.5;

F_V = repmat(F_V_pm*eye(2),[1,1,model.direction]);
% F_V = F_V_pm*eye(2);

Pol_Shape = zeros(1,model.direction);

model.birth.w = (1/nbirths)*ones(nbirths,1);
model.birth.irregularGGIW = repmat(struct('a',120,'b',10,'m',[],'P',diag([4,4,15,15]),'v',22.5,'V',85.5*eye(2),'Shape_coefficients',Pol_Shape,...
                            'ScaleEplision',65.5*eye(2),...
                            'shape',struct('num_parameter',F_a,'inverse_scale',F_b,'m',[],...
                            'deviation_matrix',F_P,'degrees_freedom',F_v,'shape_parameter',{F_V},...
                            'EDM_coefficients',model.coefficients_dilation))...
                            ,...
                            [nbirths,1]);

% model.birth.GGIW = repmat(struct('a',120,'b',10,'m',[],'P',diag([10,10,20,20]),'v',22.5,'V',85.5*eye(2)),[nbirths,1]);
model.birth.GGIW = repmat(struct('a',150,'b',10,'m',[],'P',diag([10,10,20,20]),'v',22.5,'V',90.5*eye(2)),[nbirths,1]);

model.birth.boundary = 1;
for i = 1:nbirths
    model.birth.GGIW(i).m = xstart(:,i);

    model.birth.irregularGGIW(i).m = xstart(:,i);
    model.birth.irregularGGIW(i).shape.m = xstart(:,i);
    %赋予形状的初始傅里叶参数
    model.birth.irregularGGIW(i).Shape_coefficients(floor(model.direction/2)+1) = F_V_pm/F_v_pm;
end

model.hislik = [];

%采样参数
D = F_V_pm/F_v_pm;
% blockdist为最小弦长，该弦为协方差中最大值组成圆上单位弧度(单位弧度 = 总长度/分区数量)对应的弦长
% 弦长dist = 2sin(theta/2)*R
dist = 1.2*2*sin(2*pi/model.direction)*D;
% dist = 1;
model.range_linear = Linear_equationUdist(dist, [0;0], model.direction_angle);

model.GOSPAtitles = {
    'GOSPA',
    'LocationError',
    'shapeError',
    'MissedError',
    'FalseError',
};

%% ks-pmbm
%衍生目标模型
model.deri.w = 0.001;
model.deri.GGIW = struct('a',120,'b',10,'m',[],'P',diag([4 4 20 20]),'v',22.2,'V',85.5*eye(2));
model.deri.w_death = 1;
model.ada.noise = [5 0;0 5];% 拟合
model.deri.gamma= chi2inv(Pg,model.measmodel.d);

% 剔除面积量测比过大的成分
model.falseMeasT = 0.01;
