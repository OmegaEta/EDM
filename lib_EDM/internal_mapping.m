% Calculate the displacement of internal Measurements.
function scale = internal_mapping(theta, r, boundary_points, scales)
    % Convert to polar coordinates.
    [theta_b, r_b] = cart2pol(boundary_points(:,1), boundary_points(:,2));
    
    % Construct the RBF kernel matrix.
    gamma = 0.5;
    distances = arrayfun(@(t) angular_distance(t, theta), theta_b);
    weights = exp(-gamma * distances.^2);
    
    % Normalize the weights.
    weights = weights / sum(weights);
    
    % Calculate the interpolation ratio.
    scale = dot(weights,scales);
end