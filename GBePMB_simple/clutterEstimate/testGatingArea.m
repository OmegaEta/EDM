
% figure(4);
% clf(4);
% hold on;
% 
% S = [9 0;0 4];
% center = [0;0];
% 
% range_x = [-200 200];
% range_y = [-200 200];
% 
% Z = [];
% nz = 60;
% zx = rand(nz,1)*(range_x(2)-range_x(1)) + range_x(1);
% zy = rand(nz,1)*(range_y(2)-range_y(1)) + range_y(1);
% Z = [Z [zx zy]'];
% plot(Z(1,:),Z(2,:),"r.");
% 
% 
% d = 2;
% 
% in_gate = false(size(Z,2),1);
% 
% %Take the extent into account when doing ellipsoidal gating
% S = (S + S')/2;
% 
% nu = Z - repmat(center,[1,size(Z,2)]);
% dist= sum((inv(chol(S))'*nu).^2);
% 
% %Returns measurement indices inside the gate
% gamma_size = chi2inv(0.999,d);
% in_gate(dist<gamma_size) = true;
% 
% plot(Z(1,in_gate),Z(2,in_gate),"b.");
% 
% rate = sum(in_gate)/nz;
% Area1 = rate*(range_x(2)-range_x(1))*(range_y(2)-range_y(1));
% Area2 = pi*(det(S)^0.5)*gamma_size;

Nc = [];
Nc_in_t = [];
Nc_out_t = [];
N=100;
for i = 1:N
    fprintf("%d ",i);
    [nc,nc_in_t,nc_out_t] = estimateClutter();
    Nc = [Nc;nc];
    Nc_in_t = [Nc_in_t;nc_in_t];
    Nc_out_t = [Nc_out_t;nc_out_t];
end
figure(6);
clf(6);
hold on;

plot(1:N,200*ones(N,1),'k');
plot(1:N,Nc,'r');
% plot(1:100,Nc_in_t,'y');
plot(1:N,Nc_out_t,'b');
fprintf("\nnc=%f;nc_out=%d;nc_in=%f",mean(Nc),mean(Nc_out_t),mean(Nc_out_t));

lambda_C = [];
alpha = 800;
beta = 5;
for t = 1:N
    alpha = alpha + Nc(t);
    beta = beta + 1;
    lambda_C = [lambda_C alpha/beta];

end
plot(1:N,lambda_C,"r","LineWidth",2);

%%

function [nc_hat,nc_in_hat,nc_out_hat] = estimateClutter()
% figure(5);
% clf(5);
% grid on;
% axis([-200,200,-200,200]);
% hold on;

%% 模型参数初始化
motionmodel = paraOfMotionmodel(1);
birthmodel = paraOfBirthmodel(1);
measmodel = paraOfMeasmodel(1);
%% 产生真实轨迹
T=100;
[gt,gt_t] = generateGroundtruth(1,birthmodel,motionmodel);
nt = size(gt,2);
%% 产生量测
[Z,clutter] = generateMeas(gt,measmodel,motionmodel);


I = 50;
Tc = 1;
lambda_X = [];
psg = [];
ng_in = [];
for t = 0:Tc-1
    ZZ = Z{I-t};
    GT_t = gt_t(I-t);
    % C = clutter{I};
    % 
    plot(ZZ(1,:),ZZ(2,:),"r.");
    % plot(C(1,:),C(2,:),"ko");
    IN_GATE = false(size(ZZ,2),1);
    V = [];
    for i = 1:size(GT_t.x_lambda,2)
         in_gate = gating(GT_t.xr(1:2,i),ZZ,GT_t.X(:,:,i));
         IN_GATE = IN_GATE | in_gate;
         Pg_Size = chi2inv(0.999,2);
         S = 2*GT_t.X(:,:,i);
         S =(S'+S)/2;
         v = ellipsoidalVolume(S,Pg_Size,measmodel);
         V = [V v];
    end
    lambda_X_t = sum(GT_t.x_lambda);
    ng_in_t = sum(IN_GATE);
    ng_out_t = sum(~IN_GATE);
    psg_t = sum(V)*measmodel.c_pdf;
    w_in_t = psg_t;
    w_out_t = 1-psg_t;
    
    lambda_X = [lambda_X lambda_X_t];
    psg = [psg psg_t];
    ng_in = [ng_in ng_in_t];

    if t==0
        nc_out_hat = ng_out_t/(1-psg_t);
        w_in = w_in_t;
        w_out = w_out_t;
    end
    % plot(ZZ(1,IN_GATE),ZZ(2,IN_GATE),"b.");
end

% ng_in = [ng_in_t];
% psg = [psg_t];
% lambda_X = [lambda_X_t];
nc_in_hat = optimizeClutter(ng_in,psg,lambda_X,Tc);
nc_hat = nc_out_hat*w_out+nc_in_hat*w_in;
if isnan(nc_hat)
    pause(1000);
    fprintf("error");
end

% fprintf("nc=%f,nc_in=%f,nc_out=%f\n",nc,nc_in,nc_out_t);
end
%% 
function in_gate = gating(center,Z,S)
    d = 2;
    in_gate = false(size(Z,2),1);

    %Take the extent into account when doing ellipsoidal gating
    S = (S + S')/2;
    
    nu = Z - repmat(center,[1,size(Z,2)]);
    dist= sum((inv(chol(S))'*nu).^2);
    
    %Returns measurement indices inside the gate
    gamma_size = chi2inv(0.9999,d);
    in_gate(dist<gamma_size) = true;
    
    % plot(Z(1,in_gate),Z(2,in_gate),"b.");
    
end

%% 函数
function [Z,clutter] = generateMeas(gt,measmodel,motionmodel)
% 产生量测

T = 100;
nt = size(gt,2);
Z = cell(T,1);
%target-generated measurements
for i = 1:nt
    for t = gt(i).x_bt:gt(i).x_dt
        %generate measurements only if the target is detected
        %no misdetection at time 1
        if rand < measmodel.Pd || t == 1
            nz = poissrnd(gt(i).x_lambda(t-gt(i).x_bt+1));
            Z{t} = [Z{t} mvnrnd(gt(i).xr(1:2,t-gt(i).x_bt+1)',...
                measmodel.Ch*gt(i).X(:,:,t-gt(i).x_bt+1)+measmodel.Cv,nz)'];
        end
    end
end

%append clutter
clutter = cell(T,1);
for i = 1:T
    nz = poissrnd(measmodel.c_lambda);
    zx = rand(nz,1)*(motionmodel.range_x(2)-motionmodel.range_x(1)) + motionmodel.range_x(1);
    zy = rand(nz,1)*(motionmodel.range_y(2)-motionmodel.range_y(1)) + motionmodel.range_y(1);
    clutter{i} = [zx zy]';
    Z{i} = [Z{i} [zx zy]'];
end
end

function measmodel = paraOfMeasmodel(scenario)
% 初始化量测模型
% @scenario - 场景选择
% @return - 量测模型

if scenario == 1
    %detection probability
    Pd = 0.999;
    %Poisson clutter rate
    c_lambda = 200;
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
    Ch = diag([1/4 1/4]);
    %covariance of the measurement noise
    Cv = diag([1/8 1/8]);
    
    %struct representation
    measmodel.dz = dz;
    measmodel.H = H;
    measmodel.Ch = Ch;
    measmodel.Cv= Cv;
    measmodel.Pd = Pd;
    measmodel.c_intensity = c_intensity;
    measmodel.c_lambda = c_lambda;
    measmodel.c_pdf = c_pdf;
end
end

function motionmodel = paraOfMotionmodel(scenario)
% 初始化运动模型。
% @scenario - 运动模型选择
% @return - 运动模型

if scenario==1
    %survival probability
    Ps = 0.99;
    
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
    q = 0.01;
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
end
end
function birthmodel = paraOfBirthmodel(scenario)
% 初始化出生模型。
% @scenario - 场景选择
% @return - 出生模型

if scenario==1
    %Poisson birth model with Gaussian mixture intensity
    %number of Gaussian components
    b_nG = 4;
    %Poisson intensity
    birthmodel = repmat(struct('w',log(0.05),'xr',[],'Cr',diag([5 5 4 4])^2,...
        'V',[],'v',100,'alpha',1000,'beta',100),b_nG,1);
    %specify kinematic state (x-position,y-position,x-velocity,y_velocity)
    birthmodel(1).xr = [50 50 0 0]';
    birthmodel(2).xr = [-50 50 0 0]';
    birthmodel(3).xr = [50 -50 0 0]';
    birthmodel(4).xr = [-50 50 0 0]';
    %specify extent state (orientation,two axis lengths)
    birthmodel(1).V = (500+6)*diag([1 1])^2;
    birthmodel(2).V = (500+6)*diag([1 1])^2;
    birthmodel(3).V = (500+6)*diag([1 1])^2;
    birthmodel(4).V = (500+6)*diag([1 1])^2;
    %Poisson birth rate
    b_lambda = sum(exp([birthmodel.w]));

end
end
function [gt,gt_t] = generateGroundtruth(scenario,birthmodel,motionmodel)
% 产生真实目标的运动轨迹
%total time steps
T = 100;
if scenario==1
    % 在t=1时，四个目标同时出生，并已知匀速运动到t=100；随后，在t=20，t=40，t=60均有若干个目标出生。

    %create memory to store ground truth
    %birth time, death time, target state
    gt = struct('x_bt',[],'x_dt',[],'x_lambda',[],'xr',[],'X',[]);
    
    %number of targets
    nt = 0;

    range_x=motionmodel.range_x;
    range_y=motionmodel.range_y;

    for i = 1:T-1
        %sample the number of newborn targets
        b_lambda =sum(exp([birthmodel.w]));
        nb = poissrnd(b_lambda);
        %for the first time step, make sure that at least one target is born
        if i == 1 && nb == 0
            nb = 1;
        end
        for j = 1:nb
            %number of targets increases by 1
            nt = nt + 1;
            %sample a Gaussian component in the Poisson birth intensity
            b_idx = find(rand < cumsum([0 exp([birthmodel.w])]/b_lambda),1) - 1;
            %sample an initial target state
            gt(nt).x_bt = i;
            gt(nt).x_dt = i;
            %assume fixed Poisson rate
            gt(nt).x_lambda = gamrnd(birthmodel(b_idx).alpha,1/birthmodel(b_idx).beta);
            gt(nt).xr = mvnrnd(birthmodel(b_idx).xr',birthmodel(b_idx).Cr)';
            gt(nt).X = iwishrnd(birthmodel(b_idx).V,birthmodel(b_idx).v-3);
        end
        %generate target trajectory for all the newborn targets
        for j = 1:nt
            %termintate the trajectory if the target dies or moves out of the
            %surveillance area
            %also assumes that no target dies if there is only one
            if (rand < motionmodel.Ps || nt==1) && gt(j).xr(1,end) >= range_x(1) ...
                    && gt(j).xr(1,end) <= range_x(2) ...
                    && gt(j).xr(2,end) >= range_y(1) ...
                    && gt(j).xr(2,end) <= range_y(2) && gt(j).x_dt == i
                gt(j).x_dt = i+1;
                gt(j).x_lambda = [gt(j).x_lambda gt(j).x_lambda(end)];
                %add motion noise when generating trajectory
                gt(j).xr = [gt(j).xr mvnrnd((motionmodel.Ar*gt(j).xr(:,end))',motionmodel.Cwr)'];
                gt(j).X = cat(3,gt(j).X,gt(j).X(:,:,end));
            end
        end
    end

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

end



