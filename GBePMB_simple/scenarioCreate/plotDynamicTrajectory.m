function plotDynamicTrajectory(gt_t,gt,motionmodel,figid,pt)
% 绘制目标轨迹，并按照时刻帧绘制目标扩展
% @gt - 目标集合
% @figid - 绘图隔窗编号
% @return

figure(figid);
clf(figid);
grid on
axis([-200,200,-200,200]);
hold on;
nt = size(gt,2);
T = 100;
for i = 1:nt
    %plot trajectory
    gt_ks_plot = plot(gt(i).xr(1,:),gt(i).xr(2,:),'b','linewidth',0.3);
end
% % % 
% 按时刻绘制目标扩展
for t = 1:T
    fprintf('%d ',t);
    H = [];
    %figure(2);
    for i = 1:size(gt_t(t).x_lambda,2)
        gt_es_plot = plot_extent_iw(gt_t(t).xr(1:2,i),gt_t(t).X(:,:,i),'-','r',1);
        H = [H gt_es_plot];
    end
    pause(pt);
    delete(H);
end
xlim(motionmodel.range_x)
ylim(motionmodel.range_y)
xlabel('x (m)')
ylabel('y (m)')
legend([gt_ks_plot], {'Target trajectory'});
end

