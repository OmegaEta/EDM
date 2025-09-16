function [lc_upd,bern_updated] = BGGppp_upd(ppp,z,measmodel,paras)
%PPP_UPD performs the measurement update of a PPP and computes the
%measurement likelihood

Nu = length(ppp);
l_c = zeros(Nu,1);
pd_hat = zeros(Nu,1);
for i = 1:Nu
    [l_c(i),states(i)] = ggiw(ppp(i),z,measmodel,paras);
end

for i = 1:Nu
    states(i).s = ppp(i).s + 1;
    states(i).t = ppp(i).t;
    % pd_hat(i) = log(states(i).s/(states(i).s+states(i).t));
end

if paras.Pd_estimator == 1
    for i = 1:Nu
        pd_hat(i) = states(i).s/(states(i).s+states(i).t);
    end
elseif paras.Pd_estimator == 2
    pd_hat = paras.Pd_hat;
else
    pd_hat = paras.Pd_hat;
end

[w,temp] = normalizeLogWeights([ppp.w]'+l_c+log(pd_hat));%+pd_hat
[temp,lc_upd] = normalizeLogWeights([temp size(z,2)*log(paras.clutter.c_intensity_hat)]);

%soft assignment    
[bern_updated.xr,bern_updated.Cr] = kinematic_merge(states,w);
[bern_updated.V,bern_updated.v] = extent_merge(states,w);
[bern_updated.alpha,bern_updated.beta] = gamma_merge(states,w);

[~,wid] = max(w);
labelset = [states.label];
labellastset = [states.labellast];
bern_updated.label = labelset(wid);
bern_updated.labellast = labellastset(wid);
if paras.Pd_est
    % [bern_updated.s,bern_updated.t] = beta_merge1(states,w);


    if paras.Pd_estimator == 1
        [bern_updated.s,bern_updated.t] = beta_merge_nobias(states,w);
    elseif paras.Pd_estimator == 2
        bern_updated.s = paras.Pd_hat*10;
        bern_updated.t = 10-bern_updated.s;
    else
        bern_updated.s = paras.Pd_hat*10;
        bern_updated.t = 10-bern_updated.s;
    end

else
    bern_updated.s = paras.Pd_hat*10;
    bern_updated.t = (1-paras.Pd_hat)*10;
end
%hard assignment
% [~,idx] = max(w);
% bern_updated = states(idx);
bern_updated.r = exp(temp(1));

bern_updated.isDetected = true;

global Label_id;
Label_id = Label_id+1;

if bern_updated.label>0
    % l_lastlabel = length(num2str(floor(bern_updated.label)));
    % decimal_part = floor(bern_updated.label)/(10^l_lastlabel);
    % bern_updated.label = Label_id+decimal_part;

    bern_updated.labellast = bern_updated.label;
    bern_updated.label = Label_id;
else
    
    bern_updated.label = Label_id;
    bern_updated.labellast = -1;
end

end