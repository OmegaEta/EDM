function [GGIW] = predictGGIW(GGIW,model)
GGIW_.a = GGIW(end).a/model.eta;
GGIW_.b = GGIW(end).b/model.eta;

GGIW_.m = model.motionmodel.f(GGIW(end).m);
F = model.motionmodel.F(GGIW(end).m);
GGIW_.P = F*GGIW(end).P*F' + model.motionmodel.Q;

d = 2;
M = eye(2);
e = exp(-model.Ts/model.tao);
GGIW_.v = 2*d + 2 + e*(GGIW(end).v - 2*d - 2);
GGIW_.V = e*M*GGIW(end).V*M';
%Prediction of the scaled ellipse. Ss
GGIW_.ScaleEplision = e*M*GGIW(end).ScaleEplision*M';

F_a = zeros(1,model.direction);
F_b = ones(1,model.direction);
F_P = diag([0,0,0,0]);%shape_p
F_v = zeros(1,model.direction);
F_V_pm = 0;

F_V = repmat(F_V_pm*eye(2),[1,1,model.direction]);

Pol_Shape = zeros(1,model.direction);

GGIW_ = struct(...
    'a', GGIW_.a,'b', GGIW_.b, 'm', GGIW_.m,'P', GGIW_.P,'v', GGIW_.v,'V', GGIW_.V,'ScaleEplision',GGIW_.ScaleEplision, ...
    'Shape_coefficients',Pol_Shape,...
    'shape',struct('num_parameter',F_a,'inverse_scale',F_b,'m',[],...
                            'deviation_matrix',F_P,'degrees_freedom',F_v,'shape_parameter',{F_V},...
                            'EDM_coefficients',model.coefficients_dilation));

% Shape Predict
for i = 1:model.direction
    GGIW_.shape.shape_parameter(:,:,i) = e*M*GGIW(end).shape.shape_parameter(:,:,i)*M';
end
GGIW_.shape.degrees_freedom = e*(GGIW(end).shape.degrees_freedom-2*d-2)+2*d+2;
GGIW_.shape.m = model.motionmodel.f(GGIW(end).shape.m);
GGIW_.shape.deviation_matrix =  F*GGIW(end).shape.deviation_matrix*F'+model.motionmodel.Q;
GGIW_.shape.num_parameter = GGIW(end).shape.num_parameter;
GGIW_.shape.inverse_scale = GGIW(end).shape.inverse_scale;
GGIW_.Shape_coefficients = GGIW(end).Shape_coefficients;

% % %% Some functions related to the rotating shape
% % If the length of GGIW is greater than 3
% if length(GGIW) > 3
%     % Calculate the angle theta
%     theta = 2 * (atan2(GGIW(end).m(4), GGIW(end).m(3)) - atan2(GGIW(end - 1).m(4), GGIW(end - 1).m(3)));
% % elseif length(GGIW) > 1 && length(GGIW) < 3
% %     theta = 2 * atan2(GGIW(end).m(4), GGIW(end).m(3)) - atan2(GGIW(end - 1).m(4), GGIW(end - 1).m(3));
% else
%     % If the length of GGIW does not meet the above conditions, set theta to 0
%     theta = 0;
% end
% 
% % Sample coefficient
% K = 30;
% % Generate a linearly spaced vector of angles
% Angle_detail = linspace(0, 2 * pi, model.direction * K);
% % Evaluate the Fourier series of the dilation coefficients
% dilation_coefficients = evaluateFourierSeries(GGIW(end).shape.dilation_coefficients.a, GGIW(end).shape.dilation_coefficients.b, Angle_detail);
% % Cyclically shift the data (shift k sample points to the left)
% k = round(theta / (2 * pi) * model.direction * K); % Calculate the number of shifted samples
% f_shifted = circshift(dilation_coefficients, -k);
% % Precision conversion
% [~, closest_indices] = min(abs(Angle_detail' - model.direction_angle), [], 1);
% 
% % Calculate the shifted ACn
% shifted = FFT2FS_new(f_shifted(closest_indices));
% a_shifted = shifted.a;
% b_shifted = shifted.b;
% 
% % Predict ACn using the rotation matrix
% [a_pred, b_pred] = predictShiftedACn(GGIW(end).shape.dilation_coefficients.a, GGIW(end).shape.dilation_coefficients.b, theta);
% 
% % Predict the dilation function of the shape function at the next moment based on the velocity angle
% GGIW_.shape.dilation_coefficients.a = a_pred;
% GGIW_.shape.dilation_coefficients.b = b_pred;

GGIW_.shape.EDM_coefficients.a = GGIW(end).shape.EDM_coefficients.a;
GGIW_.shape.EDM_coefficients.b = GGIW(end).shape.EDM_coefficients.b;

% 状态更新
GGIW = [GGIW;GGIW_];

end

