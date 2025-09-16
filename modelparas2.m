% load('Scenario_closeSpaced_Irregular.mat')
load('Scenario\Scenario_InterSect_Irregular3.mat');
load('Scenario\Scenario_Intersect_TruthIrregular3.mat')
scla = 1.5;
model.Scenario_range = [-240,240,-240,240]*scla;
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

%re-construct the data structure

% % generate tracks (ground truth)
% X = cell(K,1);
% E = cell(K,1);
% N = zeros(K,1);
% groundTruth = cell(K,1);
% % generate tracks (ground truth)
% for targetnum = 1:length(targetTracks)
%     for k = targetTracks(targetnum).birthTime:targetTracks(targetnum).deathTime
%         targetstate = targetTracks(targetnum).x(1:model.motionmodel.d,k-targetTracks(targetnum).birthTime+1);
%         targetextent = targetTracks(targetnum).X(:,:,k-targetTracks(targetnum).birthTime+1);
%         X{k} = [X{k} targetstate];
%         E{k} = cat(3,E{k},targetextent);
%         N(k) = N(k) + 1;
%     end
% end
% 
% for k = 1:K
%     groundTruth{k}.x = X{k};
%     groundTruth{k}.X = E{k};
% end

%target existence probability
model.Ps = 0.99;
%target detection probability
model.Pd = Scenario.detection_prob;

%range of the surveillance area
range_c = [-1 1;-1 1]*200;
%Poisson false alarm (clutter) rate
% lambda_c = Scenario.false_alarm_rate;
lambda_c = 140;
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
nbirths = 10;
model.nbirths = nbirths;
xstart = zeros(model.motionmodel.d,nbirths);

% xstart(:,1) = [-55  245 0 0];%1
% xstart(:,2) = [   0 245 0 0];%2
% xstart(:,3) = [ 80  245 0 0];%3
% xstart(:,4) = [ 245  125 0 0];%4
% 
% xstart(:,5) = [ -135  -245 0 0];%7
% xstart(:,6) = [   -0  -245 0 0];%8
% xstart(:,7) = [  40  -245 0 0];%9
% xstart(:,8) = [-245  -85 0 0];%10
% xstart(:,9) = [-245   -0 0 0];%11
% 
% xstart(:,10) = [-247   27  0 0];%12

xstart(:,1) = [-55  245 0 0]*scla;%1
xstart(:,2) = [   0 245 0 0]*scla;%2
xstart(:,3) = [ 80  245 0 0]*scla;%3
xstart(:,4) = [ 245  125 0 0]*scla;%4

xstart(:,5) = [ -135  -245 0 0]*scla;%7
xstart(:,6) = [   -0  -245 0 0]*scla;%8
xstart(:,7) = [  40  -245 0 0]*scla;%9
xstart(:,8) = [-245  -45 0 0]*scla;%10
xstart(:,9) = [-245   -0 0 0]*scla;%11

xstart(:,10) = [-245   25  0 0]*scla;%12

%Birth model
d = 2;
model.direction = 36;
model.direction_angle = linspace(0, 2*pi, model.direction);
%膨胀参数
model.coefficients_dilation.a = [1 zeros(1,model.direction)];model.coefficients_dilation.a(end) = [];%膨胀函数
model.coefficients_dilation.b = zeros(1,model.direction+1);model.coefficients_dilation.b(end) = [];

F_a = 8; F_a = repmat(F_a,1,model.direction);
F_b = 6; F_b = repmat(F_b,1,model.direction);
F_P = diag([4,4,4,4]);%shape_p
F_v_pm = 22.5; F_v = repmat(F_v_pm,1,model.direction);
F_V_pm = 56.5;

F_V = repmat(F_V_pm*eye(2),[1,1,model.direction]);
% F_V = F_V_pm*eye(2);

Pol_Shape = zeros(1,model.direction);

model.birth.w = (1/nbirths)*ones(nbirths,1);
model.birth.irregularGGIW = repmat(struct('a',120,'b',10,'m',[],'P',diag([16,16,12,12]),'v',20.5,'V',90.5*eye(2),'Shape_coefficients',Pol_Shape,...
                            'ScaleEplision',65.5*eye(2),...
                            'shape',struct('num_parameter',F_a,'inverse_scale',F_b,'m',[],...
                            'deviation_matrix',F_P,'degrees_freedom',F_v,'shape_parameter',{F_V},...
                            'EDM_coefficients',model.coefficients_dilation))...
                            ,...
                            [nbirths,1]);

% model.birth.GGIW = repmat(struct('a',120,'b',10,'m',[],'P',diag([10,10,20,20]),'v',22.5,'V',80.5*eye(2)),[nbirths,1]);
model.birth.GGIW = repmat(struct('a',120,'b',10,'m',[],'P',diag([16,16,12,12]),'v',22.5,'V',90.5*eye(2)),[nbirths,1]);

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
dist = 1*2*sin(pi/model.direction)*D;
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
model.deri.GGIW = struct('a',120,'b',10,'m',[],'P',diag([16 16 20 20]),'v',22.2,'V',90.5*eye(2));
model.deri.w_death = 1;
model.ada.noise = [5 0;0 5];% 拟合
model.deri.gamma= chi2inv(Pg,model.measmodel.d);

% 剔除面积量测比过大的成分
model.falseMeasT = 0.01;
