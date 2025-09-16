function [angle] = DirectionTransFunc(t,t_b,i,a,s)
% 方向过渡函数 27个目标 
% t时刻 t_b初始时刻 i方向模式  a初始朝向角度  s尺度

t1  = (t-t_b)/s;
[k,am] = DirectionBasicFunc(t1,i);
angle = atand(k);
angle = angle+a+am;

end

