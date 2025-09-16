function [lc_upd,bern_updated] = BGGppp_upd1(ppp,z,measmodel,paras)
%PPP_UPD performs the measurement update of a PPP and computes the
%measurement likelihood

Nu = length(ppp);
l_c = zeros(Nu,1);
for i = 1:Nu
    [l_c(i),states(i)] = ggiw(ppp(i),z,measmodel);
end

[w,temp] = normalizeLogWeights([ppp.w]'+l_c);
[temp,lc_upd] = normalizeLogWeights([temp size(z,2)*log(paras.clutter.c_intensity_hat)]);

%soft assignment    
[bern_updated.xr,bern_updated.Cr] = kinematic_merge(states,w);
[bern_updated.V,bern_updated.v] = extent_merge(states,w);
[bern_updated.alpha,bern_updated.beta] = gamma_merge(states,w);

%hard assignment
% [~,idx] = max(w);
% bern_updated = states(idx);
bern_updated.r = exp(temp(1));
bern_updated.isDetected = true;
end