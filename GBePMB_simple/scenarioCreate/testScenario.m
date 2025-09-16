clc;clear;
rng('default');
%% 模型参数初始化
motionmodel = paraOfMotionmodel(1);
birthmodel = paraOfBirthmodel(1);
measmodel = paraOfMeasmodel(1);
%% 产生真实轨迹
T=100;
[gt,gt_t] = generateGroundtruth(1,birthmodel,motionmodel);
nt = size(gt,2);
%% 产生量测
Z = generateMeas(gt,measmodel,motionmodel);
%% 绘制真实轨迹
plot_enable = false;
if plot_enable
    plotTrajectory(gt,motionmodel,1);
end

if plot_enable
    plotDynamicTrajectory(gt_t,gt,motionmodel,2);
end
%% 绘制量测
if plot_enable
    plotDynamicMeas(Z,3);
end 