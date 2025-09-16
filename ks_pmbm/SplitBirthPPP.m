function [Bern,lik] = SplitBirthPPP(meas_cell,log_W_deri,model,time)
%
if isempty(meas_cell) % ÁãÁ¿²â¸üÐÂ
    Bern.r = 0; 
    Bern.GGIW = struct('a',0,'b',1,'m',[1000;1000;0;0],'v',10,'P',diag([1 1 1 1]),'V',diag([1 1])); 
    lik = 0; 
    Bern.t_birth = 0;
    Bern.t_death = 0;
    Bern.w_death = 0.9;   
    return;
end

w = model.deri.w;
[ggiw] = adaMeas(meas_cell,model);
[Bern,lik] = detectionPPP(w,ggiw,meas_cell',model);
Bern.GGIW.P = diag([25,25,25,25]);
N = size(meas_cell,1);
lik_deri = N*log_W_deri;
lik = lik + lik_deri;
Bern.r = model.deri.intensity;
Bern.t_birth = time;
Bern.t_death = time;
Bern.w_death = model.deri.w_death;
lik = lik + log(Bern.r) + log(model.Pd) + log(Bern.w_death(end));
end

