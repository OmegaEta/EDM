 function measmodel = paraOfMeasmodel(scenario)
% 初始化量测模型
% @scenario - 场景选择
% @return - 量测模型

if scenario == 1
    %detection probability
    Pd = 0.9;
    %Poisson clutter rate
    c_lambda = 80;
    %surveillance area
    range_x = [-300 300];
    range_y = [-300 300];
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
    Ch = diag([1/2 1/2]); % Ch * X + Cv 为量测分布的协方差，其中X为扩展矩阵。
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
elseif scenario==2
    %detection probability
    Pd = 0.95;
    %Poisson clutter rate
    c_lambda = 140;
    %surveillance area
    range_x = [-500 500];
    range_y = [-500 500];
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
end
end

