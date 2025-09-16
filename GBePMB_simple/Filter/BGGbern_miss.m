function [l_missed,bern_missed] = BGGbern_miss(bern,measmodel,paras)
%BERN_MISS performs the misdetection update of a Bernoulli and computes the
%misdetection likelihood

% Ps = 0.99;
if paras.Pd_estimator == 1
    Pd = bern.s/(bern.s+bern.t);
elseif paras.Pd_estimator == 2
    Pd = paras.Pd_hat;
else
    Pd = paras.Pd_hat;
end


% Pd = measmodel.Pd;

%the sensor fails to the detect the target
Qd_1 = 1 - Pd;
%the sensor detects the target, but the target generates 0 measurement
Qd_2 = Pd*(bern.beta/(bern.beta+1))^bern.alpha;

%effective misdetection probability
Qd = Qd_1 + Qd_2;

temp = 1-bern.r+bern.r*Qd;
l_missed = log(temp);

bern_missed = bern;
bern_missed.r = bern.r*Qd/temp;
bern_missed.beta = 1/(Qd_1/Qd/bern.beta+Qd_2/Qd/(bern.beta+1));
% if bern.isD
bern_missed.isDetected = false;

if paras.Pd_est
    % bern_missed.s = bern.s;
    % bern_missed.t = bern.t+1;
    % r = bern.r;

    r =1;
    a = bern.alpha;
    b = bern.beta;
    temp1 = (1-r) + r*(1-Pd) + r*Pd*(b/(b+1))^a;
    bern_missed.s = bern.s*(r*(1-Pd))/temp1+(bern.s+1)*((1-r)+r*Pd*(b/(b+1))^a)/temp1;
    bern_missed.t = (bern.t+1)*(r*(1-Pd))/temp1+bern.t*((1-r) +r*Pd*(b/(b+1))^a)/temp1;
    
    % bern_missed.s = bern.s;
    % bern_missed.t = bern.t + 1;




else
    bern_missed.s = bern.s;
    bern_missed.t = bern.t;
end
bern_missed.label = bern.label;
bern_missed.labellast = bern.labellast;
end


