function ppp_missed = BGGppp_miss(ppp,measmodel,paras)
%PPP_MISS performs the misdetection update of a PPP

%the sensor fails to the detect the target
Qd_1 = 1 - measmodel.Pd;

ppp_missed = ppp;
for i = 1:length(ppp_missed)
    if paras.Pd_estimator == 1
        Pd = ppp(i).s/(ppp(i).s+ppp(i).t);
    elseif paras.Pd_estimator == 2
        Pd = paras.Pd_hat;
    else
        Pd = paras.Pd_hat;
    end
    % Pd = paras.Pd_hat;

    %the sensor detects the target, but the target generates 0 measurement
    Qd_2 = Pd*(ppp(i).beta/(ppp(i).beta+1))^ppp(i).alpha;
    Qd = Qd_1 + Qd_2;
    ppp_missed(i).w = ppp(i).w + log(Qd);
    ppp_missed(i).beta = 1/(Qd_1/Qd/ppp(i).beta+Qd_2/Qd/(ppp(i).beta+1));
    
    if paras.Pd_est
        % ppp_missed(i).s = ppp(i).s*Qd_1/Qd+(ppp(i).s+1)*Qd_2/Qd;
        % ppp_missed(i).t = (ppp(i).t+1)*Qd_1/Qd+ppp(i).t*Qd_2/Qd;


        if paras.Pd_estimator == 1
            ppp_missed(i).s = ppp(i).s;
            ppp_missed(i).t = ppp(i).t + 1;
        elseif paras.Pd_estimator == 2
            ppp_missed(i).s = ppp(i).s;
            ppp_missed(i).t = ppp(i).t;
        else
            ppp_missed(i).s = ppp(i).s;
            ppp_missed(i).t = ppp(i).t;
        end


    else
        ppp_missed(i).s = ppp(i).s;
        ppp_missed(i).t = ppp(i).t;

    end

    ppp_missed(i).isDetected = false;
    ppp_missed(i).label = ppp(i).label;
    ppp_missed(i).labellast = ppp(i).labellast;
end



end

