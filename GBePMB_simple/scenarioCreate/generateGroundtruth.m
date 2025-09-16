function [gt,gt_t,card] = generateGroundtruth(scenario,birthmodel,motionmodel,T)
% 产生真实目标的运动轨迹

if scenario==1
    %create memory to store ground truth
    %birth time, death time, target state
    gt = struct('x_bt',[],'x_dt',[],'x_lambda',[],'xr',[],'X',[]);
    
    %number of targets
    nt = 0;

    range_x=motionmodel.range_x;
    range_y=motionmodel.range_y;

    for i = 1:T-1
        %sample the number of newborn targets
        b_lambda =sum(exp([birthmodel.w]));
        nb = poissrnd(b_lambda);
        %for the first time step, make sure that at least one target is born
        if i == 1 && nb == 0
            nb = 1;
        end
        for j = 1:nb
            %number of targets increases by 1
            nt = nt + 1;
            %sample a Gaussian component in the Poisson birth intensity
            b_idx = find(rand < cumsum([0 exp([birthmodel.w])]/b_lambda),1) - 1;
            %sample an initial target state
            gt(nt).x_bt = i;
            gt(nt).x_dt = i;
            %assume fixed Poisson rate
            gt(nt).x_lambda = gamrnd(birthmodel(b_idx).alpha,1/birthmodel(b_idx).beta);
            gt(nt).xr = mvnrnd(birthmodel(b_idx).xr',birthmodel(b_idx).Cr)';
            gt(nt).X = iwishrnd(birthmodel(b_idx).V,birthmodel(b_idx).v-3);
        end
        %generate target trajectory for all the newborn targets
        for j = 1:nt
            %termintate the trajectory if the target dies or moves out of the
            %surveillance area
            %also assumes that no target dies if there is only one
            randnum = rand;
            % if(randnum > motionmodel.Ps)
            %     pause(1);
            % end
            if (randnum < motionmodel.Ps || nt==1) && gt(j).xr(1,end) >= range_x(1) ...
                    && gt(j).xr(1,end) <= range_x(2) ...
                    && gt(j).xr(2,end) >= range_y(1) ...
                    && gt(j).xr(2,end) <= range_y(2) && gt(j).x_dt == i
                gt(j).x_dt = i+1;
                gt(j).x_lambda = [gt(j).x_lambda gt(j).x_lambda(end)];
                %add motion noise when generating trajectory
                gt(j).xr = [gt(j).xr mvnrnd((motionmodel.Ar*gt(j).xr(:,end))',motionmodel.Cwr)'];
                gt(j).X = cat(3,gt(j).X,gt(j).X(:,:,end));
            end
        end
    end
elseif scenario==2
      %create memory to store ground truth
    %birth time, death time, target state
    gt = struct('x_bt',[],'x_dt',[],'x_lambda',[],'xr',[],'X',[]);
    
    %number of targets
    nt = 0;

    range_x=motionmodel.range_x;
    range_y=motionmodel.range_y;

    for i = 1:T-1
        %sample the number of newborn targets
        b_lambda =sum(exp([birthmodel.w]));
        nb = poissrnd(b_lambda);
        %for the first time step, make sure that at least one target is born
        if i == 1 && nb == 0
            nb = 1;
        end
        for j = 1:nb
            %number of targets increases by 1
            nt = nt + 1;
            %sample a Gaussian component in the Poisson birth intensity
            b_idx = find(rand < cumsum([0 exp([birthmodel.w])]/b_lambda),1) - 1;
            %sample an initial target state
            gt(nt).x_bt = i;
            gt(nt).x_dt = i;
            %assume fixed Poisson rate
            gt(nt).x_lambda = gamrnd(birthmodel(b_idx).alpha,1/birthmodel(b_idx).beta);
            gt(nt).X = iwishrnd(birthmodel(b_idx).V,birthmodel(b_idx).v-3);
            if b_idx==1
                if rand<0.33
                    gt(nt).xr = mvnrnd([50 50 -2 0],diag([10 10 2 2])^2)';
                elseif rand<0.66
                    gt(nt).xr = mvnrnd([50 50 0 -2],diag([10 10 2 2])^2)';
                else
                    gt(nt).xr = mvnrnd([50 50 -2 -2],diag([10 10 2 2])^2)';
                end
            elseif b_idx==2
                if rand<0.33
                    gt(nt).xr = mvnrnd([-50 50 2 0],diag([10 10 2 2])^2)';
                elseif rand<0.66
                    gt(nt).xr = mvnrnd([-50 50 0 -2],diag([10 10 2 2])^2)';
                else
                    gt(nt).xr = mvnrnd([-50 50 2 -2],diag([10 10 2 2])^2)';
                end
            elseif b_idx==3
                if rand<0.33
                    gt(nt).xr = mvnrnd([-50 -50 2 0],diag([10 10 2 2])^2)';
                elseif rand<0.66
                    gt(nt).xr = mvnrnd([-50 -50 0 2],diag([10 10 2 2])^2)';
                else
                    gt(nt).xr = mvnrnd([-50 -50 2 2],diag([10 10 2 2])^2)';
                end
            elseif b_idx==4
                if rand<0.33
                    gt(nt).xr = mvnrnd([50 -50 -2 0],diag([10 10 2 2])^2)';
                elseif rand<0.66
                    gt(nt).xr = mvnrnd([50 -50 0 2],diag([10 10 2 2])^2)';
                else
                    gt(nt).xr = mvnrnd([50 -50 -2 2],diag([10 10 2 2])^2)';
                end
            else 
                gt(nt).xr = mvnrnd(birthmodel(b_idx).xr',birthmodel(b_idx).Cr)';
            end
            
        end
        %generate target trajectory for all the newborn targets
        for j = 1:nt
            %termintate the trajectory if the target dies or moves out of the
            %surveillance area
            %also assumes that no target dies if there is only one
            randnum = rand;
            % if(randnum > motionmodel.Ps)
            %     pause(1);
            % end
            if (randnum < motionmodel.Ps || nt==1) && gt(j).xr(1,end) >= range_x(1) ...
                    && gt(j).xr(1,end) <= range_x(2) ...
                    && gt(j).xr(2,end) >= range_y(1) ...
                    && gt(j).xr(2,end) <= range_y(2) && gt(j).x_dt == i
                gt(j).x_dt = i+1;
                gt(j).x_lambda = [gt(j).x_lambda gt(j).x_lambda(end)];
                %add motion noise when generating trajectory
                gt(j).xr = [gt(j).xr mvnrnd((motionmodel.Ar*gt(j).xr(:,end))',motionmodel.Cwr)'];
                gt(j).X = cat(3,gt(j).X,gt(j).X(:,:,end));
            end
        end
    end
elseif scenario==3

     % 存活的场景
    range_x = [-200 200];
    range_y = [-200 200];
    %time interval
    Ts = 1;
    %transition matrix for kinematic state
    Ar = [1 0 Ts 0;
    0 1 0 Ts;
    0 0 1 0;
    0 0 0 1];
    %create memory to store ground truth
    %birth time, death time, target state
    gt = struct('x_bt',[],'x_dt',[],'x_lambda',[],'xr',[],'X',[]);

    nx = 1;
    gt(nx).x_bt = 1;
    gt(nx).x_dt = T;
    x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
    gt(nx).x_lambda = 5*ones(1,x_dbt);
    gt(nx).xr = [-90;90;0;2];
    gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
    for i= 1:x_dbt-1
        x_k1 = Ar*gt(nx).xr(:,end);
        gt(nx).xr = [gt(nx).xr x_k1];
    end
elseif scenario==4
    % 存活的场景
    range_x = [-200 200];
    range_y = [-200 200];
    %time interval
    Ts = 1;
    %transition matrix for kinematic state
    Ar = [1 0 Ts 0;
    0 1 0 Ts;
    0 0 1 0;
    0 0 0 1];
    %create memory to store ground truth
    %birth time, death time, target state
    gt = struct('x_bt',[],'x_dt',[],'x_lambda',[],'xr',[],'X',[]);

    Nx = 10;
    for nx = 1:Nx
        gt(nx).x_bt = 10*(nx-1)+1;
        gt(nx).x_dt = T;
        x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
        gt(nx).x_lambda = 5*ones(1,x_dbt);
        gt(nx).xr = [-85;85;2.5;-2.5];
        gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
        for i= 1:x_dbt-1
            x_k1 = Ar*gt(nx).xr(:,end);
            gt(nx).xr = [gt(nx).xr x_k1];
        end
    end
    Nx = 20;
    for nx = 11:Nx
        gt(nx).x_bt = 10*(nx-10-1)+1;
        gt(nx).x_dt =90;
        x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
        gt(nx).x_lambda = 5*ones(1,x_dbt);
        gt(nx).xr = [-90;-90;2.5;2.5];
        gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
        for i= 1:x_dbt-1
            x_k1 = Ar*gt(nx).xr(:,end);
            gt(nx).xr = [gt(nx).xr x_k1];
        end
    end

elseif scenario==5
     % 存活的场景
    range_x = [-200 200];
    range_y = [-200 200];
    %time interval
    Ts = 1;
    %transition matrix for kinematic state
    Ar = [1 0 Ts 0;
    0 1 0 Ts;
    0 0 1 0;
    0 0 0 1];
    %create memory to store ground truth
    %birth time, death time, target state
    gt = struct('x_bt',[],'x_dt',[],'x_lambda',[],'xr',[],'X',[]);
    
    nx = 0;

    nx = nx + 1;
    gt(nx).x_bt = 1;
    gt(nx).x_dt = T;
    x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
    gt(nx).x_lambda = 5*ones(1,x_dbt);
    gt(nx).xr = [-85;85;2.5;-2.5];
    gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
    for i= 1:x_dbt-1
        x_k1 = Ar*gt(nx).xr(:,end);
        gt(nx).xr = [gt(nx).xr x_k1];
    end

    nx = nx + 1;
    gt(nx).x_bt = 1;
    gt(nx).x_dt = 85;%90
    x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
    gt(nx).x_lambda = 5*ones(1,x_dbt);
    gt(nx).xr = [-80;-80;2.5;2.5];
    gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
    for i= 1:x_dbt-1
        x_k1 = Ar*gt(nx).xr(:,end);
        gt(nx).xr = [gt(nx).xr x_k1];
    end

    nx = nx + 1;
    gt(nx).x_bt = 50;
    gt(nx).x_dt = 70;
    x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
    gt(nx).x_lambda = 5*ones(1,x_dbt);
    gt(nx).xr = [80;20;0;4.2];
    gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
    for i= 1:x_dbt-1
        x_k1 = Ar*gt(nx).xr(:,end);
        gt(nx).xr = [gt(nx).xr x_k1];
    end

    nx = nx + 1;
    gt(nx).x_bt = 20;
    gt(nx).x_dt = 100;%90 100 50
    x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
    gt(nx).x_lambda = 5*ones(1,x_dbt);
    gt(nx).xr = [-120;-120;3;1];
    gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
    for i= 1:x_dbt-1
        x_k1 = Ar*gt(nx).xr(:,end);
        gt(nx).xr = [gt(nx).xr x_k1];
    end

elseif scenario==6
         % 存活的场景
    range_x = [-200 200];
    range_y = [-200 200];
    %time interval
    Ts = 1;
    %transition matrix for kinematic state
    Ar = [1 0 Ts 0;
    0 1 0 Ts;
    0 0 1 0;
    0 0 0 1];
    %create memory to store ground truth
    %birth time, death time, target state
    gt = struct('x_bt',[],'x_dt',[],'x_lambda',[],'xr',[],'X',[]);
    
    nx = 0;

    nx = nx + 1;
    gt(nx).x_bt = 1;
    gt(nx).x_dt = T;
    x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
    gt(nx).x_lambda = 5*ones(1,x_dbt);
    gt(nx).xr = [-110;100;2.5;-2.5];
    gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
    for i= 1:x_dbt-1
        x_k1 = Ar*gt(nx).xr(:,end);
        gt(nx).xr = [gt(nx).xr x_k1];
    end

    nx = nx + 1;
    gt(nx).x_bt = 20;
    gt(nx).x_dt = 85;
    x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
    gt(nx).x_lambda = 5*ones(1,x_dbt);
    gt(nx).xr = [-85;-80;3.5;3.5];
    gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
    for i= 1:x_dbt-1
        x_k1 = Ar*gt(nx).xr(:,end);
        gt(nx).xr = [gt(nx).xr x_k1];
    end

    nx = nx + 1;
    gt(nx).x_bt = 20;
    gt(nx).x_dt = 100;
    x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
    gt(nx).x_lambda = 5*ones(1,x_dbt);
    gt(nx).xr = [-120;-120;3;1];
    gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
    for i= 1:x_dbt-1
        x_k1 = Ar*gt(nx).xr(:,end);
        gt(nx).xr = [gt(nx).xr x_k1];
    end
    
    nx = nx + 1;
    gt(nx).x_bt = 50;
    gt(nx).x_dt = 75;
    x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
    gt(nx).x_lambda = 5*ones(1,x_dbt);
    gt(nx).xr = [80;20;0;3.5];
    gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
    for i= 1:x_dbt-1
        x_k1 = Ar*gt(nx).xr(:,end);
        gt(nx).xr = [gt(nx).xr x_k1];
    end
elseif scenario==7
         % 存活的场景
    range_x = [-200 200];
    range_y = [-200 200];
    %time interval
    Ts = 1;
    %transition matrix for kinematic state
    Ar = [1 0 Ts 0;
    0 1 0 Ts;
    0 0 1 0;
    0 0 0 1];
    %create memory to store ground truth
    %birth time, death time, target state
    gt = struct('x_bt',[],'x_dt',[],'x_lambda',[],'xr',[],'X',[]);
    
    nx = 0;

    nx = nx + 1;
    gt(nx).x_bt = 20;
    gt(nx).x_dt = 60;
    x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
    gt(nx).x_lambda = 5*ones(1,x_dbt);
    gt(nx).xr = [-100;100;2.5;-2.5];
    gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
    for i= 1:x_dbt-1
        x_k1 = Ar*gt(nx).xr(:,end);
        gt(nx).xr = [gt(nx).xr x_k1];
    end

    nx = nx + 1;
    gt(nx).x_bt = 68;
    gt(nx).x_dt = 100;
    x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
    gt(nx).x_lambda = 5*ones(1,x_dbt);
    gt(nx).xr = [20;-20;2.5;-2.5];
    gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
    for i= 1:x_dbt-1
        x_k1 = Ar*gt(nx).xr(:,end);
        gt(nx).xr = [gt(nx).xr x_k1];
    end


    % nx = nx + 1;
    % gt(nx).x_bt = 30;
    % gt(nx).x_dt = 100;
    % x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
    % gt(nx).x_lambda = 5*ones(1,x_dbt);
    % gt(nx).xr = [-100;100;2.5;-2.5];
    % gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
    % for i= 1:x_dbt-1
    %     x_k1 = Ar*gt(nx).xr(:,end);
    %     gt(nx).xr = [gt(nx).xr x_k1];
    % end
    
    
elseif scenario==8
         % 存活的场景
    range_x = [-200 200];
    range_y = [-200 200];
    %time interval
    Ts = 1;
    %transition matrix for kinematic state
    Ar = [1 0 Ts 0;
    0 1 0 Ts;
    0 0 1 0;
    0 0 0 1];
    %create memory to store ground truth
    %birth time, death time, target state
    gt = struct('x_bt',[],'x_dt',[],'x_lambda',[],'xr',[],'X',[]);
    
    nx = 0;

    nx = nx + 1;
    gt(nx).x_bt = 20;
    gt(nx).x_dt = 60;
    x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
    gt(nx).x_lambda = 5*ones(1,x_dbt);
    gt(nx).xr = [-100;100;2.5;-2.5];
    gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
    for i= 1:x_dbt-1
        x_k1 = Ar*gt(nx).xr(:,end);
        gt(nx).xr = [gt(nx).xr x_k1];
    end

    nx = nx + 1;
    gt(nx).x_bt = 68;
    gt(nx).x_dt = 100;
    x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
    gt(nx).x_lambda = 5*ones(1,x_dbt);
    gt(nx).xr = [20;-20;2.5;-2.5];
    gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
    for i= 1:x_dbt-1
        x_k1 = Ar*gt(nx).xr(:,end);
        gt(nx).xr = [gt(nx).xr x_k1];
    end

    nx = nx + 1;
    gt(nx).x_bt = 1;
    gt(nx).x_dt = 100;
    x_dbt = gt(nx).x_dt-gt(nx).x_bt+1;
    gt(nx).x_lambda = 5*ones(1,x_dbt);
    gt(nx).xr = [-100;-100;2.5;2.5];
    gt(nx).X = repmat(((200+6)*diag([1 1])^2)/106,1,1,x_dbt);
    for i= 1:x_dbt-1
        x_k1 = Ar*gt(nx).xr(:,end);
        gt(nx).xr = [gt(nx).xr x_k1];
    end
elseif scenario==9
 %create memory to store ground truth
    %birth time, death time, target state
    gt = struct('x_bt',[],'x_dt',[],'x_lambda',[],'xr',[],'X',[]);
    
    %number of targets
    nt = 0;

    range_x=motionmodel.range_x;
    range_y=motionmodel.range_y;

    i = 1;
    %sample the number of newborn targets
    b_lambda = 8;
    nb = poissrnd(b_lambda);
    %for the first time step, make sure that at least one target is born
    if i == 1 && nb == 0
        nb = 1;
    end
    for j = 1:nb
        %number of targets increases by 1
        nt = nt + 1;
        %sample a Gaussian component in the Poisson birth intensity
        b_idx = find(rand < cumsum([0 exp([birthmodel.w])]/sum(exp([birthmodel.w]))),1) - 1;
        %sample an initial target state
        gt(nt).x_bt = i;
        gt(nt).x_dt = i;
        %assume fixed Poisson rate
        gt(nt).x_lambda = gamrnd(birthmodel(b_idx).alpha,1/birthmodel(b_idx).beta);
        gt(nt).xr = mvnrnd(birthmodel(b_idx).xr',birthmodel(b_idx).Cr)';
        gt(nt).X = iwishrnd(birthmodel(b_idx).V,birthmodel(b_idx).v-3);
    end
    %generate target trajectory for all the newborn targets

    i = 50;
    %sample the number of newborn targets
    b_lambda = 4;
    nb = poissrnd(b_lambda);
    %for the first time step, make sure that at least one target is born
    if i == 1 && nb == 0
        nb = 1;
    end
    for j = 1:nb
        %number of targets increases by 1
        nt = nt + 1;
        %sample a Gaussian component in the Poisson birth intensity
        b_idx = find(rand < cumsum([0 exp([birthmodel.w])]/sum(exp([birthmodel.w]))),1) - 1;
        %sample an initial target state
        gt(nt).x_bt = i;
        gt(nt).x_dt = i;
        %assume fixed Poisson rate
        gt(nt).x_lambda = gamrnd(birthmodel(b_idx).alpha,1/birthmodel(b_idx).beta);
        gt(nt).xr = mvnrnd(birthmodel(b_idx).xr',birthmodel(b_idx).Cr)';
        gt(nt).X = iwishrnd(birthmodel(b_idx).V,birthmodel(b_idx).v-3);
    end
    %generate target trajectory for all the newborn targets



    for i = 1:99
        for j = 1:nt
            %termintate the trajectory if the target dies or moves out of the
            %surveillance area
            %also assumes that no target dies if there is only one
            randnum = rand;
            % if(randnum > motionmodel.Ps)
            %     pause(1);
            % end
            if (randnum < motionmodel.Ps || nt==1) && gt(j).xr(1,end) >= range_x(1) ...
                    && gt(j).xr(1,end) <= range_x(2) ...
                    && gt(j).xr(2,end) >= range_y(1) ...
                    && gt(j).xr(2,end) <= range_y(2) && gt(j).x_dt == i
                gt(j).x_dt = i+1;
                gt(j).x_lambda = [gt(j).x_lambda gt(j).x_lambda(end)];
                %add motion noise when generating trajectory
                gt(j).xr = [gt(j).xr mvnrnd((motionmodel.Ar*gt(j).xr(:,end))',motionmodel.Cwr)'];
                gt(j).X = cat(3,gt(j).X,gt(j).X(:,:,end));
            end
        end
    end


end




gt_t = repmat(struct('x_lambda',[],'xr',[],'X',[]),T,1);
nt = size(gt,2);
for t = 1:T
    for i = 1:nt
        if gt(i).x_bt<=t && gt(i).x_dt>=t
            gt_t(t).x_lambda = [gt_t(t).x_lambda gt(i).x_lambda(t-gt(i).x_bt+1)];
            gt_t(t).xr = [gt_t(t).xr gt(i).xr(:,t-gt(i).x_bt+1)];
            gt_t(t).X =  cat(3,gt_t(t).X,gt(i).X(:,:,t-gt(i).x_bt+1));
        end
    end
end


%cardinality of multi-target states
card = zeros(T,1);
for i = 1:T
    for j = 1:nt
        if gt(j).x_bt <= i && gt(j).x_dt >= i
            card(i) = card(i) + 1;
        end
    end
end
end

