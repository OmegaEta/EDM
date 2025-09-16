function ppp_predicted = BGGppp_pred(ppp,motionmodel,birthmodel,paras)
%PPP_PRED performs the prediction of a PPP

ppp_predicted = ppp;

for i = 1:length(ppp)
    ppp_predicted(i).w = ppp(i).w + log(motionmodel.Ps);
    
    ppp_predicted(i).xr = motionmodel.Ar*ppp(i).xr;
    ppp_predicted(i).Cr = motionmodel.Ar*ppp(i).Cr*motionmodel.Ar'+motionmodel.Cwr;
    
    ppp_predicted(i).v = 6 + exp(-motionmodel.Ts/motionmodel.tau)*(ppp(i).v-6);
    ppp_predicted(i).V = exp(-motionmodel.Ts/motionmodel.tau)*ppp(i).V;

    ppp_predicted(i).alpha = ppp(i).alpha/motionmodel.eta;
    ppp_predicted(i).beta = ppp(i).beta/motionmodel.eta;
    
    if paras.Pd_est
        ro = 1.01;
        mu = ppp(i).s/(ppp(i).s+ppp(i).t);
        s = ro*ppp(i).s*ppp(i).t/((ppp(i).s+ppp(i).t)^2)/(ppp(i).s+ppp(i).t+1);
        ppp_predicted(i).s = (mu*(1-mu)/s-1)*mu;
        ppp_predicted(i).t = (mu*(1-mu)/s-1)*(1-mu);
    else
        ppp_predicted(i).s =  ppp(i).s;
        ppp_predicted(i).t =  ppp(i).t;
    end

    % ppp_predicted(i).s = ppp(i).s;
    % ppp_predicted(i).t = ppp(i).t;
end

ppp_predicted = [ppp_predicted;birthmodel];

end

