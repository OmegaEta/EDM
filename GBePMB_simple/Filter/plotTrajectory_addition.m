function plotTrajectory_addition(gt,motionmodel,figid)
% 或者轨迹，并每五个时刻绘制一次目标扩展
% @gt - 目标集合
% @figid - 绘图隔窗编号
% @return
figure(figid);
clf(figid);
grid on;
box on;
axis([-200,200,-200,200]);
hold on;
nt = size(gt,2);
typeset =['r','g','b','c'];
start_pos = zeros(2,nt);
end_pos =  zeros(2,nt);
trajectory_handle = [];

%plot trajectory
i=1;
TH= plot(gt(i).xr(1,:),gt(i).xr(2,:),"r",'linewidth',1,'DisplayName',strcat('target',num2str(i)));
trajectory_handle = [trajectory_handle TH];
%plot extent for every 5 time steps
for j = 1:12:(gt(i).x_dt - gt(i).x_bt + 1)
    gt_es_plot = plot_extent_iw(gt(i).xr(1:2,j),gt(i).X(:,:,j),'-','k',1);
end
start_pos(:,i) = gt(i).xr(1:2,1);
end_pos(:,i) = gt(i).xr(1:2,end);

i=2;
TH= plot(gt(i).xr(1,:),gt(i).xr(2,:),"r",'linewidth',1,'DisplayName',strcat('target',num2str(i)));
%plot extent for every 5 time steps
for j = 1:12:(gt(i).x_dt - gt(i).x_bt + 1)
    gt_es_plot = plot_extent_iw(gt(i).xr(1:2,j),gt(i).X(:,:,j),'-','k',1);
end
start_pos(:,i) = gt(i).xr(1:2,1);
end_pos(:,i) = gt(i).xr(1:2,end);

i=3;
TH= plot(gt(i).xr(1,:),gt(i).xr(2,:),"b",'linewidth',1,'DisplayName',strcat('target',num2str(i)));
trajectory_handle = [trajectory_handle TH];
%plot extent for every 5 time steps
for j = 1:12:(gt(i).x_dt - gt(i).x_bt + 1)
    gt_es_plot = plot_extent_iw(gt(i).xr(1:2,j),gt(i).X(:,:,j),'-','k',1);
end
start_pos(:,i) = gt(i).xr(1:2,1);
end_pos(:,i) = gt(i).xr(1:2,end);


start_handle = plot(start_pos(1,:),start_pos(2,:),strcat('m','^'),'MarkerSize',10,'LineWidth',1);
end_handle = plot(end_pos(1,:),end_pos(2,:),strcat('m',"diamond"),'MarkerSize',10,'LineWidth',1);

xlim(motionmodel.range_x)
ylim(motionmodel.range_y)
xlabel('x (m)')
ylabel('y (m)')
legend([trajectory_handle gt_es_plot start_handle end_handle], ...
    {'Target1','Target2','Target extent','Birth position','Death position'});
end

