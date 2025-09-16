function [Bern,lik] = detectionPPPforUpdate(wp,GGIWp,C,model)

%Perform the initial EDM, assuming there is initially EDM
C_ = internal_expansion(C, model.coefficients_dilation);

[GGIW_c,lik_c] = updateGGIWforPPP(GGIWp,C_',model);
 
%& Initial shape estimation to Sc
for i = 1:length(GGIWp)
    [shape_coefficients,newshape] = estimate_shapeforPPP(GGIWp(i,1),C,model);
    GGIW_c(i,1).Shape_coefficients = shape_coefficients;%Update the initial shape function.
    GGIW_c(i,1).shape = newshape;%Update the initial shape boundary.
end
%&

w_c = lik_c + wp + log(model.Pd);
[w_hat,Bern.GGIW] = GGIW_merge_wrap(w_c,GGIW_c);

[~,Bern.GGIW.shape] = GGIW_shape_merge_wrap(w_c,[GGIW_c.shape]);

num_meas = size(C,2);

if num_meas > 1
    Bern.r = 1;
    lik = log(w_hat);
else
    Bern.r = w_hat/(w_hat+model.lambda_fa);
    lik = log(w_hat+model.lambda_fa);
end



end

