function birthmodel = paraOfBirthmodel(scenario)
% 初始化出生模型。
% @scenario - 场景选择
% @return - 出生模型

scla=1.5;
if scenario==1
    %Poisson birth model with Gaussian mixture intensity
    %number of Gaussian components
    b_nG = 4;
    %Poisson intensity
    birthmodel = repmat(struct('w',log(0.25),'xr',[],'Cr',diag([10 10 2 2])^2,...
        'V',[],'v',100,'alpha',1000,'beta',100,'isDetected',false,'s',7,'t',3,'label',-1),b_nG,1);
    %specify kinematic state (x-position,y-position,x-velocity,y_velocity)
    birthmodel(1).xr = [50 50 0 0]';
    birthmodel(2).xr = [-50 50 0 0]';
    birthmodel(3).xr = [-50 -50 0 0]';
    birthmodel(4).xr = [50 -50 0 0]';
    %specify extent state (orientation,two axis lengths)
    birthmodel(1).V = (100+6)*diag([1 1])^2;
    birthmodel(2).V = (100+6)*diag([1 1])^2;
    birthmodel(3).V = (100+6)*diag([1 1])^2;
    birthmodel(4).V = (100+6)*diag([1 1])^2;
    %Poisson birth rate
    b_lambda = sum(exp([birthmodel.w]));
elseif scenario==2
    %Poisson birth model with Gaussian mixture intensity
    %number of Gaussian components
    b_nG = 1;
    %Poisson intensity
    birthmodel = repmat(struct('w',log(0.25),'xr',[],'Cr',diag([4 4 2 2])^2,...
        'V',[],'v',100,'alpha',1000,'beta',100,...
        'isDetected',false,'s',6,'t',4,'label',-1,'labellast',-1),b_nG,1);
    %specify kinematic state (x-position,y-position,x-velocity,y_velocity)
    birthmodel(1).xr = [0 0 0 0]';

    %specify extent state (orientation,two axis lengths)
    birthmodel(1).V = (200+6)*diag([1 1])^2;


    %Poisson birth rate
    b_lambda = sum(exp([birthmodel.w]));
elseif scenario==3
    %Poisson birth model with Gaussian mixture intensity
    %number of Gaussian components
    b_nG = 5;
    %Poisson intensity
    birthmodel = repmat(struct('w',log(0.05),'xr',[],'Cr',diag([10 10 4 4])^2,...
        'V',[],'v',100,'alpha',1000,'beta',100,'isDetected',false),b_nG,1);
    %specify kinematic state (x-position,y-position,x-velocity,y_velocity)
    birthmodel(1).xr = [50 50 0 0]';
    birthmodel(2).xr = [-50 50 0 0]';
    birthmodel(3).xr = [-50 -50 0 0]';
    birthmodel(4).xr = [50 -50 0 0]';
    birthmodel(5).xr = [0 -0 0 0]';
    %specify extent state (orientation,two axis lengths)
    birthmodel(1).V = (150+6)*diag([1 1])^2;
    birthmodel(2).V = (150+6)*diag([1 1])^2;
    birthmodel(3).V = (150+6)*diag([1 1])^2;
    birthmodel(4).V = (150+6)*diag([1 1])^2;
    birthmodel(5).V = (100+6)*diag([1 1])^2;

    birthmodel(5).alpha = 600;
    birthmodel(5).w = log(0.03);
    %Poisson birth rate
    b_lambda = sum(exp([birthmodel.w]));
elseif scenario==4
    %Poisson birth model with Gaussian mixture intensity
    %number of Gaussian components
    b_nG = 4;
    %Poisson intensity
    birthmodel = repmat(struct('w',log(0.05),'xr',[],'Cr',diag([10 10 2 2])^2,...
        'V',[],'v',100,'alpha',1000,'beta',100,...
        'isDetected',false,'s',6,'t',4,'label',-1,'labellast',-1),b_nG,1);
    %specify kinematic state (x-position,y-position,x-velocity,y_velocity)
    birthmodel(1).xr = [50 50 0 0]';
    birthmodel(2).xr = [-50 50 0 0]';
    birthmodel(3).xr = [-50 -50 0 0]';
    birthmodel(4).xr = [50 -50 0 0]';
    %specify extent state (orientation,two axis lengths)
    birthmodel(1).V = (200+6)*diag([1 1])^2;
    birthmodel(2).V = (200+6)*diag([1 1])^2;
    birthmodel(3).V = (200+6)*diag([1 1])^2;
    birthmodel(4).V = (200+6)*diag([1 1])^2;

    %Poisson birth rate
    b_lambda = sum(exp([birthmodel.w]));
elseif scenario==5
    %Poisson birth model with Gaussian mixture intensity
    %number of Gaussian components
    b_nG = 4;
    %Poisson intensity
    birthmodel = repmat(struct('w',log(0.05),'xr',[],'Cr',diag([10 10 4 4])^2,...
        'V',[],'v',100,'alpha',1000,'beta',100,'isDetected',false),b_nG,1);
    %specify kinematic state (x-position,y-position,x-velocity,y_velocity)
    birthmodel(1).xr = [50 50 0 0]';
    birthmodel(2).xr = [-50 50 0 0]';
    birthmodel(3).xr = [-50 -50 0 0]';
    birthmodel(4).xr = [50 -50 0 0]';
    %specify extent state (orientation,two axis lengths)
    birthmodel(1).V = (150+6)*diag([1 1])^2;
    birthmodel(2).V = (150+6)*diag([1 1])^2;
    birthmodel(3).V = (150+6)*diag([1 1])^2;
    birthmodel(4).V = (150+6)*diag([1 1])^2;
    %Poisson birth rate
    b_lambda = sum(exp([birthmodel.w]));

elseif scenario==6 % 特定场景
    %Poisson birth model with Gaussian mixture intensity
    %number of Gaussian components
    b_nG = 3;
    %Poisson intensit4
    birthmodel = repmat(struct('w',log(0.01),'xr',[],'Cr',diag([4 4 15 15])^2,...
        'V',[],'v',22.5,'alpha',120,'beta',10,'isDetected',false,'label',-1,'labellast',-1),b_nG,1);
    %specify kinematic state (x-position,y-position,x-velocity,y_velocity)
    birthmodel(1).xr = [-192 75 0 0]';
    birthmodel(2).xr = [-192 10 0 0]';
    birthmodel(3).xr = [-192 -30 0 0]';
    %specify extent state (orientation,two axis lengths)
    birthmodel(1).V = (85.5)*diag([1 1])^2;
    birthmodel(2).V = (85.5)*diag([1 1])^2;
    birthmodel(3).V = (85.5)*diag([1 1])^2;
    %Poisson birth rate
    b_lambda = sum(exp([birthmodel.w]));


elseif scenario==7 % 特定场景
    %Poisson birth model with Gaussian mixture intensity
    %number of Gaussian components
    b_nG = 10;
    %Poisson intensit4
    birthmodel = repmat(struct('w',log(0.1),'xr',[],'Cr',diag([16 16 12 12])^2,...
        'V',[],'v',20.5,'alpha',220,'beta',10,'isDetected',false,'label',-1,'labellast',-1),b_nG,1);
    %specify kinematic state (x-position,y-position,x-velocity,y_velocity)
    birthmodel(1).xr = ([-55  245 0 0]*scla)';
    birthmodel(2).xr = ([   0 245 0 0]*scla)';
    birthmodel(3).xr = ([ 80  245 0 0]*scla)';
    birthmodel(4).xr = ([ 245  125 0 0]*scla)';
    birthmodel(5).xr = ([ -135  -245 0 0]*scla)';
    birthmodel(6).xr = ([   -0  -245 0 0]*scla)';
    birthmodel(7).xr = ([  40  -245 0 0]*scla)';
    birthmodel(8).xr = ([-245  -45 0 0]*scla)';
    birthmodel(9).xr = ([-245   -0 0 0]*scla)';
    birthmodel(10).xr = ([-245   25  0 0]*scla)';
    %specify extent state (orientation,two axis lengths)
    birthmodel(1).V = (90.5)*diag([1 1])^2;
    birthmodel(2).V = (90.5)*diag([1 1])^2;
    birthmodel(3).V = (90.5)*diag([1 1])^2;
    birthmodel(4).V = (90.5)*diag([1 1])^2;
    birthmodel(5).V = (90.5)*diag([1 1])^2;
    birthmodel(6).V = (90.5)*diag([1 1])^2;
    birthmodel(7).V = (90.5)*diag([1 1])^2;
    birthmodel(8).V = (90.5)*diag([1 1])^2;
    birthmodel(9).V = (90.5)*diag([1 1])^2;
    birthmodel(10).V = (90.5)*diag([1 1])^2;
    %Poisson birth rate
    b_lambda = sum(exp([birthmodel.w]));

elseif scenario==8 % 特定场景
    %Poisson birth model with Gaussian mixture intensity
    %number of Gaussian components
    b_nG = 1;
    %Poisson intensit4
    birthmodel = repmat(struct('w',log(0.01),'xr',[],'Cr',diag([5 5 2 2])^2,...
        'V',[],'v',100,'alpha',500,'beta',100,'isDetected',false),b_nG,1);
    %specify kinematic state (x-position,y-position,x-velocity,y_velocity)
    birthmodel(1).xr = [-100 100 0 0]';

    %specify extent state (orientation,two axis lengths)
    birthmodel(1).V = (200+6)*diag([1 1])^2;


    %Poisson birth rate
    b_lambda = sum(exp([birthmodel.w]));
end
end

