clc;clear;

figure(2);
clf(2);
hold on;
states = repmat(struct('alpha',1000,'beta',100),2,1);
states(1).alpha = 900;
states(1).beta = 110;
ws = [0.5 0.5];
[alpha,beta] = gamma_merge(states,ws);

alpha_wsum = [1000 900]*ws';
beta_wsum = [100 110]*ws';
pd_wsum = makedist('Gamma','a',alpha_wsum,'b',1/beta_wsum);
pd_wsum.plot;

pds = [];
for i = 1:length(states)
    pd = makedist('Gamma','a',states(i).alpha,'b',1/states(i).beta);
    pds = [pds pd];
    % pd.plot;
end
plotPDFwsum(pds,ws);
pd_merge = makedist('Gamma','a',alpha,'b',1/beta);
pd_merge.plot;

function plotPDFwsum(pds,ws)
x_range = 5:0.1:15;
pd_prange = zeros(1,length(x_range));
for i = 1:length(ws)
    pdi_prange = pdf(pds(i),x_range);
    pd_prange = pd_prange + ws(i).*pdi_prange;
end

plot(x_range,pd_prange,"k--","LineWidth",1.5);
end