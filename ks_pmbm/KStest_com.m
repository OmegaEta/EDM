function [TRUE] = KStest_com(XY,alpha,model)
% clarify
TRUE = 1;

d = 2;
N_XY = size(XY,1);
if N_XY<=1
    TRUE = 1;
    return;
end
% GGIW = adaMeas(XY,model); %拟合量测

XY_bar = mean(XY,1);
temp = (XY - XY_bar);
Sigma2_ = (temp'*temp)/(N_XY-1);
Sigma2_ = 0.5*(Sigma2_+Sigma2_');
mu_ = XY_bar;

%旋转量测与均值
[eig_v,Sigma] = eig(Sigma2_);
W = zeros(N_XY,2);
mu = (eig_v*mu_')';
SigmaX = sqrt(Sigma(1));
SigmaY = sqrt(Sigma(end));
for i = 1:N_XY
    W(i,:) = (eig_v*XY(i,:)')';
end

cum = 0.4;
N_cum = round(N_XY*cum);

%对x轴的KS检验
[d_cumx_max,delta_x] = KStestSingleAxis(cum,W(:,1),mu(1),SigmaX);

%对y轴KS检验
[d_cumy_max,delta_y] = KStestSingleAxis(cum,W(:,2),mu(2),SigmaY);

delta = max([delta_x delta_y]);
%%
D = KStable(N_XY,alpha);
uncertain = D*(N_cum^0.5);
if delta > D || abs(d_cumx_max)>uncertain || abs(d_cumy_max)>uncertain
    TRUE = 0;
end
if delta<D*0.4
    TRUE = 1;
end
end