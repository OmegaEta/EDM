
N = 3; %目标个数
T = model.K;
lamda = [15 15 10];
lambda_c = model.lambda_c;
%出生位置
xy0 = [-150 75;-150 -75;-12 0];

%存活时间
t = [1 100;1 100;45 90];

%步长
l = 3.6;
v3 = 0.45; %衍生目标的速度增量
v = [l l l+v3];

% 间距
D = 14; %平行位置间距

Target = cell(N,1);
Move = cell(N,1);
%实现Target和Move
for i = 1:N
    Target{i}.x = xy0(i,:)'; 
    if i<3
        Target{i}.X = [2.75 0;0 1.75];
    else
        Target{i}.X = [2.75 0;0 1.1];
    end
    Target{i}.t_birth = t(i,1);
    Target{i}.t_death = t(i,2);
    Move{i}.heading = zeros(Target{i}.t_death-Target{i}.t_birth+1,1);
    Move{i}.v = zeros(Target{i}.t_death-Target{i}.t_birth+1,1)+v(i);
end
Target{1}.x(2) = Target{1}.x(2) - D/2; 
Target{2}.x(2) = Target{2}.x(2) + D/2; 

xy = zeros(N,2);
k = zeros(N,1);
for t = 1:T
    for i = 1:N
        if t == 1 
            xy(i,:) = Target{i}.x';
        end
        if t>=Target{i}.t_birth && t<=Target{i}.t_death
            k(i) = DirectionFunc01(xy(i,1),t,i);
            Move{i}.heading(t-Target{i}.t_birth+1) = atand(k(i));
            [xy(i,1),xy(i,2)] = LocationFunc01(xy(i,1),xy(i,2),l,k(i));
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
box on;
axis([-200,200,-200,200]);
hold on;
for i = 1:N
    t_last = Target{i}.t_death-Target{i}.t_birth;
    for t = 1:t_last
%         Sigmacircle(Track{i}{t}.x(1),Track{i}{t}.x(2),Track{i}{t}.X,2,7,'--');
        Sigmacircle(Track{i}{t}.x(1),Track{i}{t}.x(2),Track{i}{t}.X,2,2);
    end
end

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
% drawMeas(Z,0,'*',4,1)
for i = 1:100
   Z{i} = Z{i}'; 
end

