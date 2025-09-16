function [Bern,lik] = detectionBern(Bern,C,model)

% % 轮廓计算
% z_bar = mean(C,2);
% for i = 1:length(GGIWp)
%     [max_dists, region_angles] = calculate_region_max_distance(C, z_bar, model.direction);
% end
% fit_original_boundary = fit_Fourier(region_angles, max_dists, model.N);
% original_boundary = evaluate_Fourier(fit_original_boundary, model.N, model.direction_angle);
% [Eplision,Round] = minEplision(original_boundary',1);


% 似然计算
[Bern.GGIW,lik] = updateGGIW(Bern.GGIW,C,model);

% %& 形状估计的更新
% [shape_coefficients,newshape] = estimate_shapeforMBM(Bern.GGIW(end),C,model);

lik = lik + log(Bern.r) + log(model.Pd) + log(Bern.w_death(end));
Bern.r = 1;

% Bern.GGIW(end).Shape_coefficients = shape_coefficients;%更新初始形状函数
% Bern.GGIW(end).shape = newshape;%更新初始形状参数

Bern.t_death = Bern.t_death(end);
Bern.w_death = 1;

end
