function [x1,y1] = LocationFunc01(x,y,l,k)
% λ�ú��� ����
angle = atand(k);

sinl = sind(angle)*l;
cosl = cosd(angle)*l;

x1 = x+cosl;
y1 = y+sinl;

end

