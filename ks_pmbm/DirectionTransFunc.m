function [angle] = DirectionTransFunc(t,t_b,i,a,s)
% ������ɺ��� 27��Ŀ�� 
% tʱ�� t_b��ʼʱ�� i����ģʽ  a��ʼ����Ƕ�  s�߶�

t1  = (t-t_b)/s;
[k,am] = DirectionBasicFunc(t1,i);
angle = atand(k);
angle = angle+a+am;

end

