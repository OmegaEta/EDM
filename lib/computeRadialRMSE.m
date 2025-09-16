function gsd = computeRadialRMSE(m1,m2,P1,P2,p1,p2)
% Radial Root-Mean-Square Difference (RRMSD) 极径均方误差

% Computes the Wasserstein Distance between to Gaussian distributions with
% mean m1 and m2 and covariance p1 and p2.

%P1是为滤波结果，P2为真实结果
% phi = zeros(length(P1),length(P2));

%N为可设定的采样精度
N = 3600;
theta = linspace(0, 2*pi, N);

r1 = evaluateFourierSeries(P1.a, P1.b, theta);
r2 = evaluateFourierSeries(P2.a, P2.b, theta);

%Computes the Radial Root-Mean-Square Difference (RRMSD) with P1 and P2
Dr = (mean((r1 - r2).^2))/(2*pi) ;
% gsd.phi = min(best_phi);
gsd.rmsd = Dr;

%different of Wasserstein Distance in direction
% D_wd = batchCovarianceDistance(p1, p2);
% gsd.gwd = norm(m1-m2,2) + real(mean(D_wd(:)))/2;
gsd.gwd = norm(m1-m2,2) + Dr;
gsd.Dr =  Dr;

end