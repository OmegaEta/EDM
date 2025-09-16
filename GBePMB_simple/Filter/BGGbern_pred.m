function bern_predicted = BGGbern_pred(bern,motionmodel,paras)
%BERN_PRED performs the prediction of a Bernoulli component

bern_predicted.r = bern.r*motionmodel.Ps;

bern_predicted.xr = motionmodel.Ar*bern.xr;
bern_predicted.Cr = motionmodel.Ar*bern.Cr*motionmodel.Ar'+motionmodel.Cwr;

bern_predicted.v = 6 + exp(-motionmodel.Ts/motionmodel.tau)*(bern.v-6);
bern_predicted.V = exp(-motionmodel.Ts/motionmodel.tau)*bern.V;

bern_predicted.alpha = bern.alpha/motionmodel.eta;
bern_predicted.beta = bern.beta/motionmodel.eta;

if paras.Pd_est
    ro = 1.01;
    mu = bern.s/(bern.s+bern.t);
    s = ro*bern.s*bern.t/((bern.s+bern.t)^2)/(bern.s+bern.t+1);
    bern_predicted.s = (mu*(1-mu)/s-1)*mu;
    bern_predicted.t = (mu*(1-mu)/s-1)*(1-mu);
else
    bern_predicted.s = bern.s;
    bern_predicted.t = bern.t;
end

bern_predicted.label = bern.label;
bern_predicted.labellast = bern.labellast;
% bern_predicted.s = bern.s;
% bern_predicted.t = bern.t;
end

