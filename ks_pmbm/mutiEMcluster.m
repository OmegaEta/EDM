function [mu,Sigma,in_gate,lik] = mutiEMcluster(Z,model)
% clarify w,mu,Sigma,k
com = 2;
d = 2;
N = 2;
H = [1 0 0 0;0 1 0 0];
adaGauss = adaMeas(Z,model); %与kstest重复计算 此处减低耦合
Sigma_GGIW = adaGauss.P;%adaGauss.V/(adaGauss.v - 2*d -2) + H*adaGauss.P*H';

if size(Sigma_GGIW,1) == 2
    sqrtP = N*sqrtm_2by2(Sigma_GGIW);
else
    sqrtP = N*sqrtm(Sigma_GGIW);
end
%%
% 在长轴两端点处放置高斯成分
phi = 0:pi/20:2*pi;
N_phi = floor(size(phi,2)/2);
xy = sqrtP*[cos(phi); sin(phi)];
xy2 = xy.*xy;
xy_sum = sum(xy2,1);
[~,idi_max1] = max(xy_sum);
if idi_max1>21 
    idi_max2 = idi_max1 - N_phi;
else
    idi_max2 = idi_max1 + N_phi;
end
components(1).m = 0.5*xy(:,idi_max1)'+adaGauss.m;
components(2).m = 0.5*xy(:,idi_max2)'+adaGauss.m;
dist = (sum((components(1).m-components(2).m).^2,2))^0.5;
components(1).Sigma = [0.5*dist 0;0 0.5*dist];
components(2).Sigma = [0.5*dist 0;0 0.5*dist];
% Sigmacircle(components(1).m(1),components(1).m(2),components(1).Sigma,2);
% Sigmacircle(components(2).m(1),components(2).m(2),components(2).Sigma,2);
%%
N = size(Z,1);
%初始化内存
mu = zeros(com,2);
Sigma = zeros(2,2,com);
w = zeros(com,1);
for i = 1:com
    mu(i,:) = components(i).m;
    Sigma(:,:,i) = components(i).Sigma;
    w(i) = 1/com;
end
GAM = zeros(com,N);
lik = zeros(N,com);
%% 
%em迭代
Lold = -inf;
for it = 1:30
    for i = 1:com
        GAM(i,:) = (w(i)*mvnpdf(Z,mu(i,:),Sigma(:,:,i)))';
    end
    GAM=GAM./repmat(sum(GAM,1),com,1);
    
    w = sum(GAM,2)/N;
    for i = 1:com
        mu(i,:) = sum(GAM(i,:)'.*Z,1)./sum(GAM(i,:),2);
        Z_mu = (Z-repmat(mu(i,:),N,1)).*(sqrt(GAM(i,:)))';
        Sigma_ = Sigma(:,:,i);
        Sigma(:,:,i) = Z_mu'*Z_mu/sum(GAM(i,:),2);
        Sigma(:,:,i) = 0.5*(Sigma(:,:,i)+Sigma(:,:,i)');
        eig_v = eig(Sigma(:,:,i));
        if any(find(eig_v<=1e-10))
            Sigma(:,:,i) = Sigma_; %如果Sigma为非正定矩阵 则保持不变
        end
    end
    
    N_ck = sum(GAM,2);
    pi_k = N_ck/N;
    for ic = 1:com
        lik(:,ic) = pi_k(ic)*mvnpdf(Z,mu(:,ic)',Sigma(:,:,ic));
    end
    Lnew= sum(log(sum(lik,2)));
    
    if i>1 && abs(Lnew-Lold)<0.01
        break
    else
        Lold = Lnew;
    end
end
lik = Lnew;
%%
in_gate = GAM(1,:)>GAM(2,:);
end

