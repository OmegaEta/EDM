function y = evaluateFourierSeries(a, b, x)
% Input:
    %   a: Fourier cosine coefficients (including a0)
    %   b: Fourier sine coefficients
    %   x: Discrete independent variable (supports vectors or scalars), x represents the angles to be calculated
% Output: s
    %   y: The value of the Fourier series at x

    % Initialize the output.
    y = a(1) * ones(size(x)); % The a0
    

    N_harmonics = min(floor(length(a)/2), length(a)-1); 
    for n = 1:N_harmonics
        if n <= length(a)-1
            an = a(n+1);
        else
            an = 0;
        end
        if n <= length(b)-1
            bn = b(n+1);
        else
            bn = 0;
        end
        y = y + an*cos(n*x) + bn*sin(n*x);
    end


end