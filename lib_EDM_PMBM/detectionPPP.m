function [Bern,lik] = detectionPPP(wp,GGIWp,C,model)

%执行初始膨胀，默认初始不膨胀
C_ = internal_expansion(C, model.coefficients_dilation);
%膨胀GGIW
[GGIW_c,lik_c] = updateGGIWforPPP(GGIWp,C_',model);

w_c = lik_c + wp + log(model.Pd);
[w_hat,Bern.GGIW] = GGIW_merge_wrapforPPP(w_c,GGIW_c);

num_meas = size(C,2);

if num_meas > 1
    Bern.r = 1;
    lik = log(w_hat);
else
    Bern.r = w_hat/(w_hat+model.lambda_fa);
    lik = log(w_hat+model.lambda_fa);
end

% lik = log(w_hat + model.lambda_fa^num_meas);


end

