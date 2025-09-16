function [a_pred, b_pred] = predictShiftedACn(a, b, alpha)
    N = length(a);
    a_pred = zeros(size(a));
    b_pred = zeros(size(b));
    
    for n = 0:N-1
        % Construct the rotation matrix.
        theta = n * alpha;
        R = [cos(theta),  sin(theta);
            -sin(theta), cos(theta)];
        
        % Apply the rotation matrix.
        orig_coeffs = [a(n+1); b(n+1)];
        new_coeffs = R * orig_coeffs;
        
        a_pred(n+1) = new_coeffs(1);
        b_pred(n+1) = new_coeffs(2);
    end
end