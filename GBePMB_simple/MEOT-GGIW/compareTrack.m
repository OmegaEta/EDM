clc;clear;
load('est1.mat');
load('est3.mat');
load('d_gospa1.mat');
load('d_gospa3.mat');
load('gt_t.mat');
load('decomposed_cost1.mat');
load('decomposed_cost3.mat');

figure(1);
clf(1);

axis([-200 200 -200 200]);
hold on;

T = 100;

for t = 1:T
    %红
    [~, ~,H1] = Sigmacircle(est1{t}.xr(1,:)',est1{t}.xr(2,:)',est1{t}.X,2,1,'-',0);
    %蓝 
    [~, ~,H3] = Sigmacircle(est3{t}.xr(1,:)',est3{t}.xr(2,:)',est3{t}.X,2,3,'--',0);
    [~, ~,H] = Sigmacircle(gt_t{t}.x(1,:)',gt_t{t}.x(2,:)',gt_t{t}.X,2,8,'-',1);
    fprintf('%d:%f %f %f %f \n',t,d_gospa1(t)-d_gospa3(t), ...
        decomposed_cost1(t).localisation-decomposed_cost3(t).localisation, ...
        decomposed_cost1(t).missed-decomposed_cost3(t).missed, ...
        decomposed_cost1(t).false-decomposed_cost3(t).false);
    
    pause(0.5);
    if (t<T)
        delete(H1);
        delete(H3);
        delete(H);
    end
end


