function X = sqrtm_2by2(A)
%SQRTM     Matrix square root.
%   X = SQRTM_2by2(A) is the principal square root of the matrix A, i.e. X*X = A.
%          
%   X is the unique square root for which every eigenvalue has nonnegative
%   real part.  If A has any eigenvalues with negative real parts then a
%   complex result is produced.  If A is singular then A may not have a
%   square root.  A warning is printed if exact singularity is detected.
%          
% Adapted for speed for 2x2 matrices from the MathWorks sqrtm.m implementation.
% Tim Bailey 2004.

%原始方法
% [Q, T] = schur(A);        % T is real/complex according to A.
% %[Q, T] = rsf2csf(Q, T);   % T is now complex Schur form.
% 
% R = zeros(2);
% 
% R(1,1) = sqrt(T(1,1));
% R(2,2) = sqrt(T(2,2));
% R(1,2) = T(1,2) / (R(1,1) + R(2,2));
% 
% X = Q*R*Q';

% 新方法
% A: 2x2 symmetric positive definite matrix
    a = A(1,1); b = A(1,2); c = A(2,2);

    t = sqrt((a + c)/2 + sqrt(((a - c)/2)^2 + b^2));
    s = sqrt((a + c)/2 - sqrt(((a - c)/2)^2 + b^2));

    % 构造正交变换矩阵
    if abs(b) > 1e-12
        v1 = [a - c; 2*b]; v1 = v1 / norm(v1);
        v2 = [-v1(2); v1(1)];
        Q = [v1, v2];
    else
        Q = eye(2);
    end

    X = Q * diag([t, s]) * Q';

end
