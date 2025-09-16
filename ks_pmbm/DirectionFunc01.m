function [k] = DirectionFunc01(x,t,i)
% 方向函数 衍生场景
k = [];
y1_ = (1/150)*x;
y2_ = -(1/150)*x;
y3_ = 0;

k = [k y1_ y2_ y3_];

k = k(i);
end

