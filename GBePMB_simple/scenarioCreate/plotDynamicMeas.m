function plotDynamicMeas(Z,figid,pt)
%
T = 100;
for t = 1:T
    t
    figure(figid);
    clf(figid);
    grid on
    axis([-200,200,-200,200]);
    hold on;
    plot_meas = plot(Z{t}(1,:),Z{t}(2,:),"*r");
    pause(pt);
    delete(plot_meas);
end
end

