function motionmodel = paraOfMotionmodel(scenario)
% 初始化运动模型。
% @scenario - 运动模型选择
% @return - 运动模型

if scenario==1
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
    range_x = [-500 500];
    range_y = [-500 500];

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
elseif scenario==2
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
    range_x = [-500 500];
    range_y = [-500 500];

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

