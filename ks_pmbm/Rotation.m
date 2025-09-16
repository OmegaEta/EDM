function [Ptrue,M]= Rotation(P,angle)
num = size(angle,2);
Ptrue = P;
for i=1:num
    M = [cos(angle(i)) -sin(angle(i));
         sin(angle(i)) cos(angle(i))];
    Ptrue(:,:,i) = M*P(:,:,i)*M';
end
    