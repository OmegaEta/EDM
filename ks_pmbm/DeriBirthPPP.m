function [PPP,gating_matrix_u,gating_matrix_u1] = DeriBirthPPP(MBM,PPP,W1,W,gating_matrix_d,gating_matrix_u1,gating_matrix_u,used_meas_d,model)
% 添加衍生ppp 和衍生gate
%%==========================新增函数=============================
N_W1 = size(W1,2);
idi_used_meas_d = find(used_meas_d);
N_W = size(W,2);
N_tracks = size(MBM.track,1);
if N_W1 ==0||N_tracks==0
    return;
end
gating_matrix_d1 = zeros(N_W1,N_tracks);
for i = 1:N_tracks
    gating_matrix_d1(:,i) = sum(gating_matrix_d{i},2);
end
% 同门限ks组（找出邻近目标和门限集）
[Bern_ks_group,meas_ks_group] = KScellPartition(N_tracks,N_W1,gating_matrix_d1);

N_ks_group = size(meas_ks_group,1);
idi_meas_ksgroup = cell(N_ks_group,1);
for i = 1:N_ks_group
    idi_meas_ksgroup{i} = find(meas_ks_group{i});
end

ksem_Z = cell(N_ks_group,1);
ksem_Z_idi = cell(N_ks_group,1);
% 将同门限ks组中的量测集，用ksem划分
for i = 1:N_ks_group
    [ksem_Z{i},~,ksem_Z_idi{i}] = KSEM_com01(W1',meas_ks_group{i},model);
end

for i = 1:N_ks_group
    I_deri = 0;
    gating_deri = [];
    used_deri = [];
    PPP_deri.w = [];
    PPP_deri.GGIW = [];
    N_Bern_ksgroup_i = size(Bern_ks_group{i},2);
    N_ksem_Z_i = size(ksem_Z{i},1);
    N_meas_simgksgroup = size(find(meas_ks_group{i}),1);
    %如果同门限ks组中，ksem划分出的量测子集数小于等于同门限ks组中的目标数
    if N_Bern_ksgroup_i>=N_ksem_Z_i
        continue;
    end
    lamda_prior = 0;
    AssoTable = -inf(N_Bern_ksgroup_i,N_ksem_Z_i);
    
    for j = 1:N_Bern_ksgroup_i
        lamda_sum = 0;
        %目标的假设数
        N_hyo = size(MBM.track{Bern_ks_group{i}(j)},1);
        %建立初始关联
        for k = 1:N_hyo
            %对于每一个同门限ks组中的每一个目标关联至不同量测子集的似然
            at = InitAssociation(MBM.track{Bern_ks_group{i}(j)}(k).Bern,ksem_Z{i},model);
            idi = AssoTable(j,:)<at;
            AssoTable(j,idi) = at(:,idi);
            lamda_sum = lamda_sum+MBM.track{Bern_ks_group{i}(j)}(k).Bern.GGIW(end).a/MBM.track{Bern_ks_group{i}(j)}(k).Bern.GGIW(end).b;
        end
        %ks组的先验平均lamda
        lamda_prior = lamda_prior + lamda_sum/N_hyo;
    end
    %拍卖算法建立初始关联
    [~,IDD] = auctionAlgorithm(AssoTable);
    IDI = logical(sum(IDD,1));
    %衍生目标索引
    idi_deri = find(sum(IDI,1)==0);
    idi_dete = find(sum(IDI,1)>0);
    %添加PPP gating
    N_deri = size(idi_deri,2);
    for j = 1:N_deri
        %ppp
        I_deri = I_deri + 1;
        XY = ksem_Z{i}{idi_deri(j)};
        
        d =2;
        in_gate = false(size(XY',2),1);
        % 同门限ks组中的量测子集被认为是衍生量测集后，还需要通过门限进一步筛选衍生量测
        for t = 1:N_Bern_ksgroup_i
            N_hyo = size(MBM.track{Bern_ks_group{i}(t)},1);
            for k = 1:N_hyo
                GGIW = MBM.track{Bern_ks_group{i}(t)}(k).Bern.GGIW(end);
                H = model.measmodel.H(GGIW.m);
                S = GGIW.V/(GGIW.v-2*d-2) + H*GGIW.P*H' + model.measmodel.R;
                S = (S + S')/2;
                nu = XY' - repmat(model.measmodel.h(GGIW.m),[1,size(XY',2)]);
                dist= sum((inv(chol(S))'*nu).^2);
                in_gate(dist<model.deri.gamma) = true;
            end
                
        end
        %没有被门限到的量测，认为是衍生量测
        XY = XY(~in_gate,:);
        N_XY = size(XY,1);
        if N_XY<3
           continue; 
        end
        XY_bar = mean(XY);
        temp = XY-XY_bar;
        Sigma2_ = (temp'*temp)/(N_XY-1);
        Sigma2_ = 0.5*(Sigma2_+Sigma2_')+[0.2 0;0 0.2];
        GGIW.a = model.deri.GGIW.a;
        GGIW.b = model.deri.GGIW.b;
        GGIW.m = [XY_bar';0;0];
        GGIW.P = diag([4 4 8 8]);
        GGIW.v = N_XY+model.deri.GGIW.v;
        GGIW.V = (Sigma2_*(N_XY-1))+model.deri.GGIW.V;
        %泊松权
        lamda_simglederi = model.deri.GGIW.a/model.deri.GGIW.b;
        p_poiss = zeros(5,1);
        for k = 0:4
            p_poiss(k+1) = poisspdf(N_meas_simgksgroup,lamda_prior+k*lamda_simglederi);
        end
        b_poiss = p_poiss/sum(p_poiss,1);
        %衍生目标数
        w_poiss = sum((0:4)'.*b_poiss,1);
        PPP_deri.w = [PPP_deri.w;log(model.deri.w * model.Ps+w_poiss)];
        PPP_deri.GGIW = [PPP_deri.GGIW;GGIW];
        %gating
        falseG = false(N_W1,1);
        falseD = false(N_W,1);
        falseG(ksem_Z_idi{i}{idi_deri(j)}) = true;
        falseD(idi_used_meas_d(falseG)) = true;
        gating_deri =  [gating_deri falseG];
        used_deri = [used_deri falseD];
    end
    %根据同门限ks组添加衍生ppp
    PPP_deri.w = [PPP_deri.w;log(model.deri.w(1))];
    GGIW = model.deri.GGIW;
    GGIW.m = [mean(W1(:,meas_ks_group{i})')';0;0]; 
    PPP_deri.GGIW = [PPP_deri.GGIW;GGIW];
    PPP.w = [PPP.w;PPP_deri.w;];
    PPP.GGIW = [PPP.GGIW;PPP_deri.GGIW];
    gating_matrix_u1 = logical([gating_matrix_u1 gating_deri sum(gating_deri,2)]);
    gating_matrix_u = logical([gating_matrix_u used_deri sum(used_deri,2)]);
    
end


end
%%=============================================================
