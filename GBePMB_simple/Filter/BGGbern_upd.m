function [l_upd,bern_updated] = BGGbern_upd(bern,z,measmodel,paras)
%BERN_UPD performs the measurement update of a Bernoulli and computes the
%measurement likelihood

[l_upd,bern_updated] = ggiw(bern,z,measmodel,paras);
bern_updated.r = 1;

if paras.Pd_estimator == 1
    Pd = bern.s/(bern.s+bern.t);
elseif paras.Pd_estimator == 2
    Pd = paras.Pd_hat;
else
    Pd = paras.Pd_hat;
end


l_upd = l_upd + log(bern.r) + log(Pd);% + Pd;
bern_updated.isDetected = true;
if paras.Pd_est
    bern_updated.s = bern.s + 1;
    bern_updated.t = bern.t;
else
    bern_updated.s = bern.s;
    bern_updated.t = bern.t;

end

bern_updated.label = bern.label;
bern_updated.labellast = bern.labellast;
end
