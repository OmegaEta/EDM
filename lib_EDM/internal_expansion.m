function [new_internal] = internal_expansion(internal, coefficients_dilation)
    %Calculate the difference.
    center = mean(internal, 2);              % 2×1
    internal_centered = internal - center;   % 2×N
    internal_centered = internal_centered';  % 转为 N×2

    % Calculate the angle and evaluate the scale.
    internal_theta = atan2(internal_centered(:,2), internal_centered(:,1));     % N×1
    internal_scale = evaluateFourierSeries(coefficients_dilation.a, coefficients_dilation.b, internal_theta); % N×1

    % Directly vectorize the scaling.
    new_internal = internal_centered .* internal_scale + center';   % N×2 .* N×1 and 1×2
end
