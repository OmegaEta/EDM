function [measureinblock_index] = irregular_PerBlockCov_gating(W,M,model)
    % Obtain the number of directions and calculate the mean.
    direction = model.direction;
    num_points = size(W, 2);
    W_ = W - M(1:2);
    
    % Initialize the index matrix.
    measureinblock_index = false(direction, num_points);
    range_linear = model.range_linear;
   
    for j = 1:direction
        % % Obtain the parameters of the current Sector.
        a = range_linear.a(:, j);
        b = range_linear.b(:, j);
        c = range_linear.c(:, j);
        cond = range_linear.plusorminus(:,:,j);
        
        % % Calculate the three inequality conditions for all points simultaneously.
        eq = [...
            a(1)*W_(1,:) + b(1)*W_(2,:) + c(1) < 0; 
            a(2)*W_(1,:) + b(2)*W_(2,:) + c(2) < 0;
            a(3)*W_(1,:) + b(3)*W_(2,:) + c(3) < 0];
        
        % Determine the points that satisfy all three conditions simultaneously.
        measureinblock_index(j,:) = all(eq == cond(:), 1);
    end
end