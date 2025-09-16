function [d_cum_max,delta] = KStestSingleAxis(cum,W,mu,Sigma)
%clarify

N_XY = size(W,1);
N_cum = round(N_XY*cum);

p_e = linspace(0,1,N_XY+1);
p_e = p_e(2:end);

p0 = normcdf(W,mu,Sigma);
p_sort = sort(p0)';
p_d = p_e - p_sort;
d_k1 = abs(p_d);
% figure(5);
% plot(sort(W),p_sort);
% hold on;
% plot(sort(W),p_e,'*');
%不确定度的累加，找出最大的不确定度区间
d_cum = zeros(N_XY-N_cum+1,1);
for i_x = 1:N_XY-N_cum+1
    for j = 1:N_cum
        d_cum(i_x) = d_cum(i_x) + p_d(i_x+j-1);
    end
end
d_cum_max = max(abs(d_cum));
% d_cumx_min = min(d_cum);

d_k2 = zeros(1,N_XY-1);
%量测过少不适合错位相减
if N_XY>=5
    d_k2 = abs(p_sort(2:end) - p_e(1:end-1));
end
delta = max([d_k1 d_k2]);
end

