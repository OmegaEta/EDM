
N = 28; %目标个数
T =100;
lamda = 15+zeros(N,1);
lambda_c = 60;%杂波数
%出生位置
xy0 =  [75 75;75 75;75 75;75 75;75 75;75 75;75 75;
        75 -75;75 -75;75 -75;75 -75;75 -75;75 -75;75 -75;
        -75 -75;-75 -75;-75 -75;-75 -75;-75 -75;-75 -75;-75 -75;
        -75 75;-75 75;-75 75;-75 75;-75 75;-75 75;-75 75];

%存活时间段
t = [1 70;
    8 100;
    15 79;
    19 50;
    25 54;
    29 100;
    35 100;
     23 100;
     2 23;
     40 90;
     45 63;
     26 70;
     22 60;
     74 86;
      6 83;
      15 86;
      24 77;
      40 100;
      44 100;
      67 94;
      76 100;
       7 100;
       16 30;
       5 25;
       38 100;
       47 100;
       59 100;
       85 99];
%速度
l = 3.6;
V = cell(N,1);
for i=1:N
    t_last = t(i,2)-t(i,1)+1;
    V{i} = l+zeros(t_last,1);
end

%方向模式
I = [1 3 3 4 5 6 7 ...
    1 2 3 4 5 6 7 ...
    1 2 3 4 5 6 7 ...
    1 2 3 4 5 6 7];
%尺度
S = [10;10;10;9;10;10;10;
    15*ones(7,1);
    10;10;10;10;10;10;10;15*ones(7,1)];
% 初始朝向
A = [-75;-73;-135;-135;45;-70;-130;
    120;210;180;230;120;-45;45;
    45;135;20;-75;-30;-20;-135;
    -90;-90;180;-20;40;220;20];
%实现Target
Target = cell(N,1);
for i=1:N
    Target{i}.x = xy0(i,:)';
    Target{i}.X = [1.5 0;0 1];
    Target{i}.t_birth = t(i,1);
    Target{i}.t_death = t(i,2);
end

%实现Move
Move = cell(N,1);
for i=1:N
    Move{i}.v = V{i};
    Move{i}.heading = zeros(Target{i}.t_death-Target{i}.t_birth+1,1);
    for t = 1:T
        if t>=Target{i}.t_birth && t<=Target{i}.t_death
            [angle] = DirectionTransFunc(t,Target{i}.t_birth,I(i),A(i),S(i));
            Move{i}.heading(t-Target{i}.t_birth+1) = angle;
        end 
    end
end


%实现Track和Z0
Z0 = cell(N,1);
Track = cell(N,1);
for i = 1:N
    Track{i} = productTrackLocation(Target{i},Move{i});
    Z0{i} = productRandomMeas(Track{i},lamda(i));
end

% 画图 groundTruth;
figure(1);
clf(1);

axis([-200,200,-200,200]);
hold on;
box on;
for i = 1:N
    t_last = Target{i}.t_death-Target{i}.t_birth;
    for t = 1:t_last
        Sigmacircle(Track{i}{t}.x(1),Track{i}{t}.x(2),Track{i}{t}.X,2,2);
        [i t]
    end
end

% for t = 1:100
%     clf(3);
%     axis([-200,200,-200,200]);
%     hold on;
%     for i = 1:N
%         if t >=Target{i}.t_birth && t<=Target{i}.t_death
%             Sigmacircle(Track{i}{t-Target{i}.t_birth+1}.x(1),Track{i}{t-Target{i}.t_birth+1}.x(2),Track{i}{t-Target{i}.t_birth+1}.X,2,1);
%         end
%     end
% %     pause(0.2);
%     [t]
% end
% t时刻 t_b初始时刻 i方向模式  a初始朝向角度  s尺度


% 杂波
Clutter = cell(100,1);
for i = 1:100
    ClutterRate = random('Poisson',lambda_c);
    Clutter{i} = unifrnd(-200,200,ClutterRate,2);
end
% 产生量测集
Z = Clutter;

%实现groundTruth
groundTruth = cell(T,1);
for i = 1:T
    groundTruth{i}.x = [];
    groundTruth{i}.X = [];
end
for i = 1:N
    t_last = Target{i}.t_death-Target{i}.t_birth;
    t_birth = Target{i}.t_birth;
    for t = 1:t_last+1
        groundTruth{t+t_birth-1}.x = [groundTruth{t+t_birth-1}.x Track{i}{t}.x];
        groundTruth{t+t_birth-1}.X = cat(3,groundTruth{t+t_birth-1}.X,Track{i}{t}.X);
        if unifrnd(0,1)>=0.05
            Z{t+t_birth-1} = [Z{t+t_birth-1};Z0{i}{t}];
        end
        
    end
end

%画图Z
% drawMeas(Z,1,'.',4,1)
for i = 1:100
   Z{i} = Z{i}'; 
end



