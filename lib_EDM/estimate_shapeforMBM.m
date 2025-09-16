% The output consists of the shape parameters and shape boundary.
% shape_coefficients-->Sc
function [shape_coefficients,newshape] = estimate_shapeforMBM(shape_r,W,model)

direction = model.direction;
z_bar = mean(W,2);
measureinblock_index = irregular_PerBlockCov(W,model);
d = 2;

shape = shape_r.shape;
measureinblock_index = logical(measureinblock_index);
irregular_card_W = sum(measureinblock_index,2);%Count the number of measurements in each direction.

epsilon = z_bar - model.measmodel.h(shape_r.m);
H = model.measmodel.H(shape_r.m);

shape_.num_parameter = [shape.num_parameter] + irregular_card_W';
shape_.inverse_scale = [shape.inverse_scale] + 1;

%Parameter initialization.
X_hat = zeros(2,2,direction);S = zeros(2,2,direction);Vs = zeros(2,2,direction);inv_sqrt_S = zeros(2,2,direction);
iS = zeros(2,2,direction);K = zeros(4,2,direction);m_ = zeros(4,1,direction);P_ = zeros(4,4,direction);
X_sqrt =  zeros(2,2,direction);S_sqrt_inv = zeros(2,2,direction);N = zeros(2,2,direction);
det_S = zeros(1,direction);

%%%%% Vectorize the structure. %%%%%
X_hat = shape.shape_parameter./reshape(shape.degrees_freedom-6,[1,1,direction]);
for i = 1:model.direction
    if irregular_card_W(i) == 0
            S(:,:,i) = H*shape.deviation_matrix*H' + model.Qd;
       else
            S(:,:,i) = X_hat(:,:,i)/irregular_card_W(i) + H*shape.deviation_matrix*H';
    end
    Vs(:,:,i) = chol(S(:,:,i));
    det_S(:,i) = prod(diag(Vs(:,:,i)))^2;
end

% % Pre-allocate.
inv_sqrt_S = zeros(2,2,direction);
iS = zeros(2,2,direction);
% Extract the components of Vs.
a = squeeze(Vs(1,1,:));  % 1×N
b = squeeze(Vs(1,2,:));
% The lower triangle in the output of `chol` should be zero, but it is retained for safety.
c = squeeze(Vs(2,1,:));
d = squeeze(Vs(2,2,:));
% Calculate the determinant.
det_Vs = a .* d - b .* c;
% Calculate the 2x2 inverse.
inv_sqrt_S(1,1,:) =  d ./ det_Vs;
inv_sqrt_S(1,2,:) = -b ./ det_Vs;
inv_sqrt_S(2,1,:) = -c ./ det_Vs;
inv_sqrt_S(2,2,:) =  a ./ det_Vs;
%  iS = inv_sqrt_S * inv_sqrt_S'
%  2×2×N, can be explicitly expanded.
iS(1,1,:) = inv_sqrt_S(1,1,:).^2 + inv_sqrt_S(1,2,:).^2;
iS(2,2,:) = inv_sqrt_S(2,1,:).^2 + inv_sqrt_S(2,2,:).^2;
iS(1,2,:) = inv_sqrt_S(1,1,:) .* inv_sqrt_S(2,1,:) + inv_sqrt_S(1,2,:) .* inv_sqrt_S(2,2,:);
iS(2,1,:) = iS(1,2,:);

for i = 1:model.direction
    K(:,:,i) = shape.deviation_matrix*H'*iS(:,:,i);
    m_(:,:,i) =shape_r.m + K(:,:,i)*epsilon;

    P_(:,:,i) = shape.deviation_matrix - K(:,:,i)*H*shape.deviation_matrix;

    X_sqrt(:,:,i) = sqrtm_2by2(X_hat(:,:,i));
    S_sqrt_inv(:,:,i) = sqrtm_2by2(iS(:,:,i));
    N(:,:,i) = X_sqrt(:,:,i)*S_sqrt_inv(:,:,i)*(epsilon*epsilon')*S_sqrt_inv(:,:,i)'*X_sqrt(:,:,i)';
end

shape_.m = mean(m_,3);
shape_.deviation_matrix = mean(P_,3);
shape_.degrees_freedom =  [shape.degrees_freedom] + irregular_card_W';

% Parameter pre-extraction.
temp = W-z_bar; %Place the measured centroid at the origin.
temp1 = reshape(temp, 2, 1, size(W, 2));  % 2×1×N
temp2 = reshape(temp, 1, 2, size(W, 2));  % 1×2×N
Z_measurements = (temp1 .* temp2);

%% Calculate the covariance of the measurement samples in each direction.
R_indirection = zeros(1,direction);  %Measurents at the polar radii in the given directions.
Z_indirection = zeros(2,2,direction);  %Measurents at the covariance in each of the given directions.
for i = 1:direction
    idx = find(measureinblock_index(i,:));  % Find the valid measurement indices in this direction.
    if ~isempty(idx)
        Z_indirection(:,:,i) = sum(Z_measurements(:,:,idx), 3);
        shape_.shape_parameter(:,:,i) = shape.shape_parameter(:,:,i) + N(:,:,i) + Z_indirection(:,:,i);

        [~, D] = eig(shape_.shape_parameter(:,:,i));
        R_indirection(i) = max(D(:)) / (shape_.degrees_freedom(i) );
    else
        shape_.shape_parameter(:,:,i) = shape.shape_parameter(:,:,i) + N(:,:,i);

        [~, D] = eig(shape_.shape_parameter(:,:,i));
        R_indirection(i) = max(D(:)) / shape_.degrees_freedom(i);
    end
end

%% 
%fit Sc
shape_coefficients = FFT2FS_new(R_indirection);
shape_.EDM_coefficients = shape_r.shape.EDM_coefficients;

newshape = shape_;

end