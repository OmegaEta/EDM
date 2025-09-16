clc;clear;
load('est2_1.mat');
load('est2_3.mat');
load('d_gospa2_1.mat');
load('d_gospa2_3.mat');
load('gt2_1.mat');
load('decomposed_cost2_1.mat');
load('decomposed_cost2_3.mat');

figure(1);
clf(1);

axis([-200 200 -200 200]);
hold on;

T = 100;

for t = 1:T
    %红
    [~, ~,H1] = Sigmacircle(est2_1{t}.xr(1,:)',est2_1{t}.xr(2,:)',est2_1{t}.X,2,1,'-',0);
    %蓝 
    [~, ~,H3] = Sigmacircle(est2_3{t}.xr(1,:)',est2_3{t}.xr(2,:)',est2_3{t}.X,2,3,'--',0);
    [~, ~,H] = Sigmacircle(gt_t2_1{t}.x(1,:)',gt_t2_1{t}.x(2,:)',gt_t2_1{t}.X,2,8,'-',1);
    fprintf('%d:%f %f %f %f \n',t,d_gospa2_1(t)-d_gospa2_3(t), ...
        decomposed_cost2_1(t).localisation-decomposed_cost2_3(t).localisation, ...
        decomposed_cost2_1(t).missed-decomposed_cost2_3(t).missed, ...
        decomposed_cost2_1(t).false-decomposed_cost2_3(t).false);
    
    pause(0.5);
    if (t<T)
        delete(H1);
        delete(H3);
        delete(H);
    end
end


