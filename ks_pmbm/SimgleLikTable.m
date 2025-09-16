function [simgleLikTable] = SimgleLikTable(berns,W,Bern_ks_group,meas_ks_group,N_berns_prior,log_W_deri,model,simgleLikTable,logidid)
% 单个量测分别关联至各个bern的ggiw似然

W1 = W(meas_ks_group,:);
Berns_ks = berns(Bern_ks_group);
N_Bernsprior_ks = size(find(Bern_ks_group<=N_berns_prior),2);

N_W1 = size(W1,1);
N_Berns_ks = size(Berns_ks,1);

if nargin == 7 %输入参数的个数
    %完全计算
    logidid = true(N_Berns_ks,1);
    simgleLikTable = zeros(N_W1,N_Berns_ks);
end
for i = 1:N_W1
    for j = 1:N_Berns_ks
        if logidid(j)
            if j<=N_Bernsprior_ks
                [~,lik] = detectionBern(Berns_ks(j),W1(i,:)',model);
                simgleLikTable(i,j) = lik;
            else
                [~,lik] = splitBirthBern(Berns_ks(j),W1(i,:),log_W_deri(j-N_Bernsprior_ks),model);
                simgleLikTable(i,j) = lik;
            end
        end
    end
end

%在出生位置更可能发生新生目标，而不是衍生目标
m = [0;0];
H = [1 0 0 0;0 1 0 0];
d = 2;
N = 0;
for i = Bern_ks_group
    if berns(i).r>0
        N = N+1;
        m = m + berns(i).GGIW(end).m(1:2);
    end
end
m = m/N;
pdf =0;
for i = 1:size(model.birth.GGIW,1)
    x = model.birth.GGIW(i).m(1:2);
    X = H*model.birth.GGIW(i).P*H' + model.birth.GGIW(i).V/(model.birth.GGIW(i).v - 2*d -2);
    pdf = pdf + mvnpdf(m,x,X);
end
if log(pdf)>-8
    simgleLikTable(:,N_Bernsprior_ks+1:end) = simgleLikTable(:,N_Bernsprior_ks+1:end)-5;
end

end

