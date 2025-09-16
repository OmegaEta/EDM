function [k,am] = DirectionBasicFunc(t,i)
% 方向基函数 27个目标

% 抛物线
if i==1
    k =1.5*(t-5);
    am = 84.2894;
elseif i==2
    k = -1.5*(t-5);
    am = -84.2894;
elseif i==3
    t = t-5;
    if t<=5
        k = -t*((25-t^2)^(-0.5));
    else
        k = -inf;
    end
    am = -90;
elseif i==4
    t = t-5;
    if t<=5
        k = t*(25-t^2)^(-0.5);
    else
        k = inf;
    end
    am = 90;
else
    k = 0;
    am = 0;
end

end

