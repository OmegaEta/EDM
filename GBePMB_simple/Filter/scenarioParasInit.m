function [paras,birthmodel,measmodel,ppp] = scenarioParasInit(scenarioParas,paras,measmodel,birthmodel,ppp)
%场景参数初始化


paras.clutter_est = scenarioParas.clutter_est;%杂波率估计器开关
paras.clutter.alpha = scenarioParas.alpha_clutter;%杂波率估计器gamma分布的alpha
paras.clutter.beta = scenarioParas.beta_clutter; %杂波率估计器gamma分布的beta
paras.clutter.lambdaC_hat = paras.clutter.alpha/paras.clutter.beta; 
paras.clutter.c_intensity_hat = paras.clutter.lambdaC_hat*measmodel.c_pdf;


paras.Pd_est = scenarioParas.Pd_est; %检测概率估计器开关
paras.Pd_estimator = 2; % 检测概率估计方法：
paras.Pd_hat = scenarioParas.Pd_hat;%检测概率初始值


measmodel.Pd =NaN;% scenarioParas.Pd_hat; %真实的检测概率
measmodel.c_lambda =NaN;% paras.clutter.alpha/paras.clutter.beta;
measmodel.c_intensity = measmodel.c_lambda*measmodel.c_pdf;



%出生模型的初始beta分布参数，用于检测概率估计器
for i=1:length(birthmodel)
    birthmodel(i).s = scenarioParas.s_Pd(i);
    birthmodel(i).t = scenarioParas.t_Pd(i);
end
for i=1:length(ppp)
    ppp(i).s = scenarioParas.s_Pd(i);
    ppp(i).t = scenarioParas.t_Pd(i);
end

end