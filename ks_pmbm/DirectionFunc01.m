function [k] = DirectionFunc01(x,t,i)
% ������ ��������
k = [];
y1_ = (1/150)*x;
y2_ = -(1/150)*x;
y3_ = 0;

k = [k y1_ y2_ y3_];

k = k(i);
end

