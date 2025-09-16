% Boundary mapping function.
function scale = boundary_mapping(theta, r_original, a, b, phi)
    %Calculate the relative rotation angle.
    theta_rel = theta - phi;
    
    %Calculate the radius of the target ellipse.
    r_target = (a*b) / sqrt((b*cos(theta_rel))^2 + (a*sin(theta_rel))^2);
    
    %Calculate the discrete Extension Difference.
    scale = r_target / r_original;
end