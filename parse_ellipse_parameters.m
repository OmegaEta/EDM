%% Ellipse parameter analysis function.
function [a, b, phi] = parse_ellipse_parameters(V, D)
    % Extract the major axis length.
%     a=sqrt(1/abs(D(1)));b=sqrt(1/abs(D(4)));
    a=sqrt(D(1));b=sqrt(D(4));

    %Calculate the rotation angle.
    phi = mod(atan2(V(2,1), V(1,1))+pi,2*pi);
end