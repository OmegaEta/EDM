% New boundary, new internal points, and the corresponding coefficients of the EDM function
function [coefficients_dilation_NEW] = shape_expansion(boundary, V, D, center)
    % Obtain the elliptical parameters from eigenvalue decomposition.
    [a, b, phi] = parse_ellipse_parameters(V, D);

    % Calculate the Differences.
    boundary_centered = boundary - center;%Contour Extension difference.
    
    % Calculate the transformation of boundary
    [theta_boundary, r_original_boundary] = cart2pol(boundary_centered(:,1), boundary_centered(:,2));
    scale_factors = arrayfun(@(t,r) boundary_mapping(t, r, a, b, phi), theta_boundary, r_original_boundary);
%     save f_orig.mat scale_factors;
    % Apply the inverse rotation transformation.
    scaled_points = scale_factors .* boundary_centered;
    new_boundary = (scaled_points')' + center;
    
    % Fit the EDM function.
    coefficients_dilation_NEW = FFT2FS_new(scale_factors);

end