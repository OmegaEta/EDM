function [TrackLocation] = productTrackLocation(Target,Move)
% 根据目标和运动模型 产生轨迹位置

t_birth = Target.t_birth;
t_death = Target.t_death;
t_last = t_death - t_birth;
TrackLocation = cell(t_last,1);
x = zeros(4,1);
x(1:2) = Target.x;
v = Move.v(1);
heading = Move.heading(1);
sindAngle=sind(heading);
cosdAngle=cosd(heading);
move=[cosdAngle*v;sindAngle*v];
x(3:4) = move;
X = Target.X;
X = Rotation(X,(heading/180)*pi);
TrackLocation{1}.x = x;
TrackLocation{1}.X = X;

for t = 2:t_last+1
    heading = Move.heading(t);
    X = Target.X;
    v = Move.v(t);
    sindAngle=sind(heading);
    cosdAngle=cosd(heading);
    move=[cosdAngle*v;sindAngle*v];
    x(1:2) = x(1:2)+move;
    x(3:4) = move;
    angled = (heading/180)*pi;
    [X,~] = Rotation(X,angled);
    TrackLocation{t}.x = x;
    TrackLocation{t}.X = X;
end

end

