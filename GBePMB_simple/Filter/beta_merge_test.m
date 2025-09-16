figure(1);
clf(1);
hold on;

% demo1
% alphas = [10 2];
% betas = [2 10];
% ws = [0.6 0.4];

% demo2
% alphas = [10 3 15];
% betas = [2 3 5];
% ws = [0.9 0.05 0.05];

%demo3
% alphas = [99 20];
% betas = [1 20 ];
% ws = [0.7 0.3];

%demo4
% alphas = [18 10];
% betas = [2 1 ];
% ws = [0.9 0.1];

%demo5
% alphas = [18 10];
% betas = [2 1 ];
% ws = [0.99 0.01];

%demo6
% alphas = [18 5];
% betas = [18 8];
% ws = [0.99 0.01];

%demo7
% alphas = [2 5];
% betas = [18 11];
% ws = [0.95 0.05];

%demo8
% alphas = [15 20];
% betas = [20 15];
% ws = [0.5 0.5];

% demo9
% alphas = [2 3 5];
% betas = [10 3 15];
% ws = [0.9 0.05 0.05];

% demo10
% alphas = [800 400];
% betas = [200 120];
% ws = [0.98 0.02];

% demo11
alphas = [50 25];
betas = [2 2];
ws = [0.98 0.02];

pds = [];
pdf_mean = [];
for i = 1:length(ws)
    pd = makedist('Beta','a',alphas(i),'b',betas(i));
    pds = [pds pd];
    pdf_mean =[pdf_mean pd.mean];
    fprintf("mean=%f\n",pd.mean);
    % pd.plot("PlotType","pdf");
end
fprintf("wmean=%f\n",ws*pdf_mean');
% multi-Beta plot
plotPDFwsum(pds,ws);


% 仅仅使用参数加权求和的Beta分布
alpha_wsum = alphas*ws';
beta_wsum = betas*ws';
pd_wsum = makedist('Beta','a',alpha_wsum,'b',beta_wsum);
% pd_wsum.plot;


% n=300;
% rs100 = randomSampling([w1,w2],[pd1 pd2],n);
% pd_fit = fitdist(rs100,'Beta');
% pd_fit.plot;
[alpha1,beta1] = beta_merge(alphas,betas,ws');

[alpha2,beta2] = beta_merge_sameMean(alphas,betas,ws');

[alpha3,beta3] = beta_merge_sameMean2(alphas,betas,ws');

figure(2);
clf(2);hold on;
plotPDFwsum(pds,ws);

pd1 = makedist('Beta','a',alpha1,'b',beta1);
pd1.plot;
fprintf("\nmean=%f\n",pd1.mean);

pd2 = makedist('Beta','a',alpha2,'b',beta2);
pd2.plot;
fprintf("\nSameMean=%f\n",pd2.mean);

pd3 = makedist('Beta','a',alpha3,'b',beta3);
pd3.plot;
fprintf("\nSameMean=%f\n",pd3.mean);
% [alpha2,beta2] = beta_merge_newtons(alphas,betas,ws);
% pd = makedist('Beta','a',alpha2,'b',beta2);
% pd.plot;

function plotPDFwsum(pds,ws)
x_range = 0:0.01:1;
pd_prange = zeros(1,length(x_range));
for i = 1:length(ws)
    pdi_prange = pdf(pds(i),x_range);
    pd_prange = pd_prange + ws(i).*pdi_prange;
end

plot(x_range,pd_prange,"k--","LineWidth",1.5);
end

function rs = randomSampling(ws,pds,n)
rs = zeros(n,1);
% n_pds = length(pds);
w_cum = cumsum(ws);
randomseed = rand(n,1);
for i = 1:n
    idx = find(randomseed(1)<w_cum,1,'first');
    rs(i) = pds(idx).random;
end

end


