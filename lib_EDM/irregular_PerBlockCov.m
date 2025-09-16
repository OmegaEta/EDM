function [measureinblock_index] = irregular_PerBlockCov(W, model)
    % Get the number of directions and calculate the mean.
    direction = model.direction;
    z_bar = mean(W, 2);
    num_points = size(W, 2);
    W_ = W - z_bar;
    
    % Initialize the index matrix
    measureinblock_index = false(direction, num_points);
    range_linear = model.range_linear;
    
    % Vectorize the calculation of partition indices for each direction.
    for j = 1:direction
        % Obtain the parameters for the current direction.
        a = range_linear.a(:, j);
        b = range_linear.b(:, j);
        c = range_linear.c(:, j);
        cond = range_linear.plusorminus(:,:,j);
        
        % Simultaneously calculate the three inequality conditions for all.
        eq = [...
            a(1)*W_(1,:) + b(1)*W_(2,:) + c(1) < 0; 
            a(2)*W_(1,:) + b(2)*W_(2,:) + c(2) < 0;
            a(3)*W_(1,:) + b(3)*W_(2,:) + c(3) < 0];
        
        % Determine the points that satisfy all three conditions simultaneously.
        measureinblock_index(j,:) = all(eq == cond(:), 1);
    end
end