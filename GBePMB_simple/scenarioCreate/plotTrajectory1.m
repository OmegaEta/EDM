function plotTrajectory1(gt,motionmodel,figid)
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
start_pos = zeros(2,nt);
end_pos =  zeros(2,nt);
for i = 1:nt
    %plot trajectory
    gt_ks_plot = plot(gt(i).xr(1,:),gt(i).xr(2,:),'b','linewidth',1);
    %plot extent for every 5 time steps
    for j = 12:12:(gt(i).x_dt - gt(i).x_bt + 1)
        gt_es_plot = plot_extent_iw(gt(i).xr(1:2,j),gt(i).X(:,:,j),'-','r',1);
    end
    start_pos(:,i) = gt(i).xr(1:2,1);
    end_pos(:,i) = gt(i).xr(1:2,end);
end

% start_handle = plot(start_pos(1,:),start_pos(2,:),strcat('m','^'),'MarkerSize',10,'LineWidth',1);
% end_handle = plot(end_pos(1,:),end_pos(2,:),strcat('m',"diamond"),'MarkerSize',10,'LineWidth',1);

xlim(motionmodel.range_x)
ylim(motionmodel.range_y)
xlabel('x (m)')
ylabel('y (m)')
legend([gt_ks_plot, gt_es_plot], {'Target trajectory','Target extent'});
end

