%% 高效率方法
%输出的是形状函数
function [shape_coefficients,newshape] = estimate_shapeforPPP(shape_r,W,model)

direction = model.direction;
z_bar = mean(W,2);
measureinblock_index = irregular_PerBlockCov(W,model);

for k = 1:length(shape_r)
    shape = shape_r(k,1).shape;
    measureinblock_index = logical(measureinblock_index); % 关键修正点
    irregular_card_W = sum(measureinblock_index,2);%求每个方向量测的个数
    
    epsilon = z_bar - model.measmodel.h(shape_r(k,1).m);
    H = model.measmodel.H(shape_r(k,1).m);
    
    shape_.num_parameter =[shape.num_parameter] + irregular_card_W';
    shape_.inverse_scale =[shape.inverse_scale] + 1;
    
    %Initialize the parameters.
    X_hat = zeros(2,2,direction);S = zeros(2,2,direction);Vs = zeros(2,2,direction);inv_sqrt_S = zeros(2,2,direction);
    iS = zeros(2,2,direction);K = zeros(4,2,direction);m_ = zeros(4,1,direction);P_ = zeros(4,4,direction);
    X_sqrt =  zeros(2,2,direction);S_sqrt_inv = zeros(2,2,direction);N = zeros(2,2,direction);
    for i = 1:model.direction
        X_hat(:,:,i) = shape.shape_parameter(:,:,i)/(shape.degrees_freedom(i)-6);
        X_hat(:,:,i) = (X_hat(:,:,i)+X_hat(:,:,i))/2;
        
        if irregular_card_W(i) == 0
                S(:,:,i) = H*shape.deviation_matrix*H';
           else
                S(:,:,i) = X_hat(:,:,i)/irregular_card_W(i) + H*shape.deviation_matrix*H';
        end
        S(:,:,i) = (S(:,:,i) + S(:,:,i)')/2;
    
        Vs(:,:,i) = chol(S(:,:,i));
        inv_sqrt_S(:,:,i) = inv(Vs(:,:,i));
        iS(:,:,i) =  inv_sqrt_S(:,:,i)*inv_sqrt_S(:,:,i)';
        
        K(:,:,i) = shape.deviation_matrix*H'*iS(:,:,i);
        m_(:,:,i) =shape_r(k,1).m + K(:,:,i)*epsilon;

        P_(:,:,i) = shape.deviation_matrix - K(:,:,i)*H*shape.deviation_matrix;

        X_sqrt(:,:,i) = sqrtm_2by2(X_hat(:,:,i));
        S_sqrt_inv(:,:,i) = sqrtm_2by2(iS(:,:,i));
        N(:,:,i) = X_sqrt(:,:,i)*S_sqrt_inv(:,:,i)*(epsilon*epsilon')*S_sqrt_inv(:,:,i)'*X_sqrt(:,:,i)';
    end
    shape_.m = mean(m_,3);
    shape_.deviation_matrix = mean(P_,3);
    shape_.degrees_freedom =  [shape.degrees_freedom]-3 + irregular_card_W';

    % Pre-extract the parameters.
    temp = W-z_bar; %% Place the centroid of the measurements at the origin.
    % %新方法
    temp1 = reshape(temp, 2, 1, size(W, 2));  %  2×1×N
    temp2 = reshape(temp, 1, 2, size(W, 2));  %  1×2×N
    Z_measurements = (temp1 .* temp2) ./ reshape(sum(measureinblock_index, 1), [1,1,size(W,2)]);
    
    %Calculate the covariance of the measurement samples in each direction.
    R_indirection = zeros(1,direction);  %% Measurements at the radial distances in the given directions.
    Z_indirection = zeros(2,2,direction);  %Measurements at the covariance in the given directions.
    for i = 1:direction
        idx = find(measureinblock_index(i,:));  % Find the valid measurement indices in this direction.
        if ~isempty(idx)
            Z_indirection(:,:,i) = sum(Z_measurements(:,:,idx), 3);
            shape_.shape_parameter(:,:,i) = shape.shape_parameter(:,:,i) + N(:,:,i) + Z_indirection(:,:,i);
    
            [~, D] = eig(shape_.shape_parameter(:,:,i));
            R_indirection(i) = max(D(:)) / shape_.degrees_freedom(i);
        else
            shape_.shape_parameter(:,:,i) = shape.shape_parameter(:,:,i) + N(:,:,i);
    
            [~, D] = eig(shape_.shape_parameter(:,:,i));
            R_indirection(i) = max(D(:)) / shape_.degrees_freedom(i);  % the expansion difference.
        end
    end
    
    % 
    %fit Sc
    coefficients_shape = FFT2FS_new(R_indirection);
    shape_coefficients = coefficients_shape;
    shape_.EDM_coefficients = shape_r.shape.EDM_coefficients;

    newshape = shape_;
end

% %画图
% theta_fit = linspace(0, 2*pi, 108);
% y = max(evaluate_Fourier(coefficients_shape.acn, theta_fit),0);
% [x_mean, y_mean] = pol2cart(theta_fit, y');
% x_mean_ = x_mean+z_bar(1);y_mean_ = y_mean+z_bar(2);
% 
% figure;
% plot(x_mean_, y_mean_, 'r-'); title('拟合的形状');hold on
% plot(W(1,:),W(2,:),'k.')

end
