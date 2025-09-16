function  V = ellipsoidalVolume(S,Pg_Size,measmodel)
%计算椭圆门限内面积
d = measmodel.dz;

if d == 1
    V = 2*(det(S)^0.5)*Pg_Size;
elseif d == 2
    V = pi*(det(S)^0.5)*Pg_Size;
else
    V = (4*pi/3)*(det(S)^0.5)*Pg_Size;
end

end

