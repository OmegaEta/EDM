
load('./MBM/ÑÜÉú³¡¾°¹ì¼£/trajectoryEstimates_improved2_3.mat');
load('./MBM/ÑÜÉú³¡¾°¹ì¼£/trajectoryEstimates_original2_3.mat');

trace_original = trajectory_original{100};
N_original = length(trace_original);

trace_improved = trajectory_improved{100};
N_improved = length(trace_improved);

figure(1)
clf(1);
axis([-200,200,-200,200]);
box on;
for i=1:N_original
    X = trace_original(i).x(1,:);
    Y = trace_original(i).x(2,:);
    line(X,Y,'Color','b');
    hold on;
end

figure(2)
clf(2);
axis([-200,200,-200,200]);
box on;
for i=1:N_improved
    X = trace_improved(i).x(1,:);
    Y = trace_improved(i).x(2,:);
    line(X,Y,'Color','r');
    hold on;
end