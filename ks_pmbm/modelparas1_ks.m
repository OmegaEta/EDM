load('Scenario\Scenario_InterSect_Irregular1.mat');
load('Scenario\Scenario_Intersect_TruthIrregular.mat')
K = 100;
model.K = K;%总时刻

% Effective window length for the gamma prediction
w_e_gamma = 15;
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
% model.Pd = Scenario.detection_prob;
model.Pd = 0.95;

%range of the surveillance area
range_c = [-1 1;-1 1]*200;
%Poisson false alarm (clutter) rate
% lambda_c = Scenario.false_alarm_rate;
lambda_c = 60;
model.lambda_c = lambda_c;
%Poisson clutter intensity
model.lambda_fa = lambda_c/prod(range_c(:,2)-range_c(:,1));

% target initial state
nbirths = 4;
xstart = zeros(model.motionmodel.d,nbirths);

%出生位置
xstart(:,1) = [-192 75 0 0];
xstart(:,2) = [-192 10 0 0];
xstart(:,3) = [-192 -30 0 0];
xstart(:,4) = [500 500 0 0];

%Birth model
d = 2;
model.birth.w = 0.25*ones(nbirths,1);
model.birth.GGIW = repmat(struct('a',1500,'b',100,'m',[],'P',diag([100,100,25,25]),'v',30,'V',72.5*eye(2)),[nbirths,1]);
for i = 1:nbirths
    model.birth.GGIW(i).m = xstart(:,i);
end
%衍生目标模型
model.deri.w = 0.001;
model.deri.GGIW = struct('a',1000,'b',100,'m',[],'P',diag([100 100 25 25]),'v',25,'V',50*eye(2));
model.deri.w_death = 1;
model.ada.noise = [5 0;0 5];% 拟合

% 剔除面积量测比过大的成分
model.falseMeasT = 0.01;

% Gating parameters
Pg = 0.94;
model.deri.gamma= chi2inv(Pg,model.measmodel.d);
Pg = 0.9994;
model.gamma= chi2inv(Pg,model.measmodel.d);
model.Qd = 1 - model.Pd*Pg;

% Thresholds
model.threshold_r = 1e-2;   %existence probability of Bernoulli component
model.threshold_u = 0.01;   %weight of mixture component in PPP 0.01
model.threshold_w = 0.005;  %weight of global hypothesis (multi-Bernoulli) 1e-2
model.threshold_s = 1e-4;   %weight of the trajectory is still alive
model.recycle = 1e-1;       %recycling threshold
model.merge = 4;            %merge threshold used to merge similar GGIWs
model.M = 100;              %cap of number of MBM components in PMBM
model.num_iterations = 3;   %controls the number of iterations used in SO
model.max_repetition = 1;   %controls the number of iterations used in SO

%extract target state from Bernoulli components with existence probability
%larger than this threshold 
model.exist_r = 0.5;        

