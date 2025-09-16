function mutiTrajectoryCompare(Z,gt_t,EST,pt,T,figid)
%

colorset = ["r","b","g"];
N_est = length(EST);
figure(figid);
clf(figid);
grid on
axis([-200,200,-200,200]);
hold on;
% 按时刻绘制目标扩展
for t = 1:T
    fprintf('%d ',t);
    H = [];
    %量测
    hZ = plot(Z{t}(1,:),Z{t}(2,:),"k.");
    H = [H hZ];
    % 真实值
    for i = 1:size(gt_t(t).x_lambda,2)
        gt_es_plot = plot_extent_iw(gt_t(t).xr(1:2,i),gt_t(t).X(:,:,i),'--',"k",1);
        H = [H gt_es_plot];
    end
    % 估计值
    for j = 1:N_est
        est = EST{j};
        est_onetime = est{t};
        for k = 1:est_onetime.card
            gt_es_plot = plot_extent_iw(est_onetime.xr(1:2,k),est_onetime.X(:,:,k),'-',colorset(j),1);
            H = [H gt_es_plot];
        end
    
    end

    pause(pt);
    delete(H);
end

end

