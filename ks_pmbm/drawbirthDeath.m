s1 =scatter([83,141],[161,161],'MarkerEdgeColor','k','Marker','^','DisplayName','死亡位置')
s1.SizeData=100
s1.LineWidth = 2
hold on
s =scatter([-150,-150],[-127.5,-72.5],'MarkerEdgeColor','k','Marker','d','DisplayName','出生位置')
s.SizeData=100
s.LineWidth = 2
legend([s,s1])
xlabel('x');ylabel('y');axis equal;box on;