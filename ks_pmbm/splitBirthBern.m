function [bern,lik] = splitBirthBern(Bern_deri,meas_cell,log_W_deri,model)
%
if Bern_deri.r==0
    bern = Bern_deri;
    lik = -inf;
    return;
end
[bern,lik] = detectionBern(Bern_deri,meas_cell',model);
N = size(meas_cell,1);
lik = lik + log_W_deri*N;

end

