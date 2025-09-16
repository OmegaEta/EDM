%产生衍生目标的场景

lambda_c = model.lambda_c;
T = model.K;
% D = 5;
d = D/2;
X = [2.75 0;0 1.75];
lamda = 15;% 2*(4*(X(1)^0.5)*(4*X(end)))^0.5;
% 目标
Target{1}.x = [-150;-75];
Target{1}.X = X;
Target{1}.t_birth = 1;
Target{1}.t_death = 17;% 17
Target{2}.x = [-150;-125];
Target{2}.X = X;
Target{2}.t_birth = 1;
Target{2}.t_death = 17; %17

Target{3}.x = [-75;-100];
Target{3}.X = X;
Target{3}.t_birth = 18;
Target{3}.t_death = 38; %21

% Target{4}.x = [25;-100];
% Target{4}.X = [1.5 0;0 1];
% Target{4}.t_birth = 39;
% Target{4}.t_death = 62; %24
% Target{5}.x = [25;-103];
% Target{5}.X = [1.5 0;0 1];
% Target{5}.t_birth = 39;
% Target{5}.t_death = 62; %24

Target{4}.x = [25;-100];
Target{4}.X = X;
Target{4}.t_birth = 39;
Target{4}.t_death = 63;% 24

Target{5}.x = [100;-25];
Target{5}.X = X;
Target{5}.t_birth = 64;
Target{5}.t_death = 83; %21

Target{6}.x = [100;75];
Target{6}.X = X;
Target{6}.t_birth = 84;
Target{6}.t_death = 100; %17
Target{7}.x = [100;75];
Target{7}.X = X;
Target{7}.t_birth = 84;
Target{7}.t_death = 100; %17

%运动轨迹
%靠近上
duration = Target{1}.t_death-Target{1}.t_birth+1;
angle = atand(3/4)/duration;
angle_cum = cumsum(zeros(duration,1)+angle)-angle;
Move{1}.heading = (360-atand(3/4)) + angle_cum;
Move{1}.v = zeros(duration,1)+(2*pi*125*atand(3/4)/360)/(duration-1); %18
%靠近下
duration = Target{2}.t_death-Target{2}.t_birth+1;
angle = atand(3/4)/duration;
angle_cum = cumsum(zeros(duration,1)+angle);
Move{2}.heading = flip(angle_cum);
Move{2}.v = zeros(duration,1)+(2*pi*125*atand(3/4)/360)/(duration-1); %18

%水平直线
Target{3}.x = Target{3}.x + [Move{1}.v(end);0];
duration = Target{3}.t_death-Target{3}.t_birth+1;
Move{3}.heading = zeros(duration,1);
Move{3}.v = zeros(duration,1)+100/(duration-1); %

% %半弧内
% angle = 90/24;
% angle_cum = cumsum(zeros(24,1)+angle)-angle;
% Move{4}.heading = angle_cum;
% Move{4}.v = 2*pi*74*0.25/24; %24
% % 半弧外
% angle = 90/24;
% angle_cum = cumsum(zeros(24,1)+angle)-angle;
% Move{5}.heading = angle_cum;
% Move{5}.v = 2*pi*75*0.25/24; %24

%半弧
Target{4}.x = Target{4}.x + [Move{1}.v(end)+Move{3}.v(end);0];
duration = Target{4}.t_death-Target{4}.t_birth+1;
angle = 90/duration;
angle_cum = cumsum(zeros(duration,1)+angle)-angle;
Move{4}.heading = angle_cum;
Move{4}.v = zeros(duration,1)+(2*pi*(75-d)*90/360)/(duration-1); %18

duration = Target{4}.t_death-Target{4}.t_birth+1;
angle = 90/duration;
angle_cum = cumsum(zeros(duration,1)+angle)-angle;
Move{5}.heading = angle_cum;
Move{5}.v = zeros(duration,1)+(2*pi*(75+d)*90/360)/(duration-1); %18

%竖直直线
Target{5}.x = Target{5}.x + [Move{1}.v(end)+Move{3}.v(end)+0.5;Move{4}.v(end)+0.75];
duration = Target{5}.t_death-Target{5}.t_birth+1;
Move{6}.heading = 90+zeros(duration,1);
Move{6}.v = zeros(duration,1)+100/(duration-1); %

%分离左
Target{6}.x = Target{6}.x + [Move{1}.v(end)+Move{3}.v(end)+0.5;Move{4}.v(end)+Move{6}.v(end)+0.75];
duration = Target{6}.t_death-Target{6}.t_birth+1;
angle = atand(3/4)/duration;
angle_cum = cumsum(zeros(duration,1)+angle)-angle;
Move{7}.heading = 90 + angle_cum;
Move{7}.v = zeros(duration,1)+(2*pi*125*atand(3/4)/360)/(duration-1); %

%分离右
Target{7}.x = Target{7}.x + [Move{1}.v(end)+Move{3}.v(end)+0.5;Move{4}.v(end)+Move{6}.v(end)+0.75];
duration = Target{7}.t_death-Target{7}.t_birth+1;
angle = atand(3/4)/duration;
angle_cum = cumsum(zeros(duration,1)+angle);
Move{8}.heading = 90 - angle_cum;
Move{8}.v = zeros(duration,1)+(2*pi*125*atand(3/4)/360)/(duration-1); %


%组合一：靠近
Track1 = productTrackLocation(Target{1},Move{1});
Track2 = productTrackLocation(Target{2},Move{2});
for i = 1:size(Track1,1) %让两个目标分开一些距离
    Track1{i}.x = Track1{i}.x + [0;d;0;0];
    Track2{i}.x = Track2{i}.x + [0;-d;0;0];
end
group1.Track{1,1} = Track1;
group1.Track{2,1} = Track2;
group1.time = [Target{1}.t_birth Target{1}.t_death;
               Target{2}.t_birth Target{2}.t_death];
    %产生量测
Z1 = productRandomMeas(Track1,lamda);
Z2 = productRandomMeas(Track2,lamda);
group1Meas.TrackMeas{1,1} = Z1;
group1Meas.TrackMeas{2,1} = Z2;
group1Meas.time = group1.time;
S.ScenarioMeas{1,1} = group1Meas;
groundTruth_.group{1,1} = group1;

%组合二：两个目标水平并行
Track1 = productTrackLocation(Target{3},Move{3});
Track2 = productTrackLocation(Target{3},Move{3});
for i = 1:size(Track1,1)
    Track1{i}.x = Track1{i}.x + [0;d;0;0];
    Track2{i}.x = Track2{i}.x + [0;-d;0;0];
end
group2.Track{1,1} = Track1;
group2.Track{2,1} = Track2;
group2.time = [Target{3}.t_birth Target{3}.t_death;
               Target{3}.t_birth Target{3}.t_death];
    %产生量测
Z1 = productRandomMeas(Track1,lamda);
Z2 = productRandomMeas(Track2,lamda);
group2Meas.TrackMeas{1,1} = Z1;
group2Meas.TrackMeas{2,1} = Z2;
group2Meas.time = group2.time;
S.ScenarioMeas{2,1} = group2Meas;
groundTruth_.group{2,1} = group2;

%情景三：两个目标并行弧
Track1 = productTrackLocation(Target{4},Move{4});
Track2 = productTrackLocation(Target{4},Move{5});
for i = 1:size(Track1,1)
    Track1{i}.x = Track1{i}.x + [0;d;0;0];
    Track2{i}.x = Track2{i}.x + [0;-d;0;0];
end
group3.Track{1,1} = Track1;
group3.Track{2,1} = Track2;
group3.time = [Target{4}.t_birth Target{4}.t_death;
               Target{4}.t_birth Target{4}.t_death];
    %产生量测
Z1 = productRandomMeas(Track1,lamda);
Z2 = productRandomMeas(Track2,lamda);
group3Meas.TrackMeas{1,1} = Z1;
group3Meas.TrackMeas{2,1} = Z2;
group3Meas.time = group3.time;
S.ScenarioMeas{3,1} = group3Meas;
groundTruth_.group{3,1} = group3;

%组合四：两个目标竖直并行
Track1 = productTrackLocation(Target{5},Move{6});
Track2 = productTrackLocation(Target{5},Move{6});
for i = 1:size(Track1,1)
    Track1{i}.x = Track1{i}.x + [-d;0;0;0];
    Track2{i}.x = Track2{i}.x + [d;0;0;0];
end
group4.Track{1,1} = Track1;
group4.Track{2,1} = Track2;
group4.time = [Target{5}.t_birth Target{5}.t_death;
               Target{5}.t_birth Target{5}.t_death];
    %产生量测
Z1 = productRandomMeas(Track1,lamda);
Z2 = productRandomMeas(Track2,lamda);
group4Meas.TrackMeas{1,1} = Z1;
group4Meas.TrackMeas{2,1} = Z2;
group4Meas.time = group4.time;
S.ScenarioMeas{4,1} = group4Meas;
groundTruth_.group{4,1} = group4;

%组合五：分离
Track1 = productTrackLocation(Target{6},Move{7});
Track2 = productTrackLocation(Target{7},Move{8});
for i = 1:size(Track1,1) %让两个目标分开一些距离
    Track1{i}.x = Track1{i}.x + [-d;0;0;0];
    Track2{i}.x = Track2{i}.x + [d;0;0;0];
end
group5.Track{1,1} = Track1;
group5.Track{2,1} = Track2;
group5.time = [Target{6}.t_birth Target{6}.t_death;
               Target{7}.t_birth Target{7}.t_death];
    %产生量测
Z1 = productRandomMeas(Track1,lamda);
Z2 = productRandomMeas(Track2,lamda);
group5Meas.TrackMeas{1,1} = Z1;
group5Meas.TrackMeas{2,1} = Z2;
group5Meas.time = group5.time ;
S.ScenarioMeas{5,1} = group5Meas;
groundTruth_.group{5,1} = group5;

% 杂波
Clutter = cell(100,1);
for i = 1:100
    ClutterRate = random('Poisson',lambda_c);
    Clutter{i} = unifrnd(-200,200,ClutterRate,2);
end


% 产生量测集
Z = Clutter;
% Z = cell(100,1);
for i = 1:size(S.ScenarioMeas,1)
    for j = 1:size(S.ScenarioMeas{i}.time,1)
        k = 0;
        for t = S.ScenarioMeas{i}.time(j,1):S.ScenarioMeas{i}.time(j,2)
            k = k+1;
            if unifrnd(0,1)>0.95
                continue;
            end
            Z{t} = [Z{t};S.ScenarioMeas{i}.TrackMeas{j}{k}];
        end
    end
end

groundTruth = cell(100,1);
for i = 1:100
    groundTruth{i}.x = [];
    groundTruth{i}.X = [];
end
figure(1);

clf(1);
box on;
axis([-200,200,-200,200]);
hold on;
for i = 1:size(groundTruth_.group,1)
    for j = 1:size(groundTruth_.group{i}.Track,1)
        time =  groundTruth_.group{i}.time(j,1):groundTruth_.group{i}.time(j,2);
        for t = 1:size(groundTruth_.group{i}.Track{j},1)
            x = groundTruth_.group{i}.Track{j}{t}.x;
            X = groundTruth_.group{i}.Track{j}{t}.X;
            %画图
             Sigmacircle(x(1),x(2),X,2,7,'--');
             Sigmacircle(x(1),x(2),X,2,2);
            groundTruth{time(t)}.x = [groundTruth{time(t)}.x x];
            groundTruth{time(t)}.X = cat(3,groundTruth{time(t)}.X,X);
        end
    end
end
    
%画图
% drawMeas(Z,0,'*',4,1)
for i = 1:100
   Z{i} = Z{i}'; 
end
