function [ksem_Z,C_t,ksem_Z_idi] = KSEM_com01(W,idi,model)
%==================��������===================

XY = W(idi,:);
if size(XY,1)==0
    ksem_Z = cell(1,1);
    ksem_Z_idi = cell(1);
    C_t = [];
    return;
end
alpha = 0.1;%������ˮƽ 
N_splitMax = 100; %��󻮷���
Max_iteration = 50;%��������

cell_spliting = cell(N_splitMax,1);%�����鼯
cell_splited = cell(N_splitMax,1);%������ɼ�
N_splited = 1;
cell_spliting{1} = XY;
C_f =[];%α�Ӳ���

% ks����
for i = 1:Max_iteration
    %cell_spliting���Ƿ���Ԫ�أ�����У��ҵ�����
    idi0 = zeros(N_splitMax,1);
    for i0 = 1:N_splitMax %�ҳ������鼯������
        if ~isempty(cell_spliting{i0})
            idi0(i0) = 1;
        end
    end
    idi = find(idi0==1);
    idi_emp = find(idi0==0);
    if isempty(idi) %����⼯Ϊ��
        break;
    end
    
    for i1 = idi'         
        [TRUE] = KStest_com(cell_spliting{i1},alpha,model);
        if TRUE == 1
            % �������
            cell_splited{N_splited} = cell_spliting{i1};
            N_splited = N_splited + 1;
            cell_spliting{i1} = [];
        else
            % ���鲻���� em���¾���
            [mu,Sigma,in_gate,lik] = mutiEMcluster(cell_spliting{i1},model);

            XY_1 = cell_spliting{i1}(in_gate==1,:);
            XY_2 = cell_spliting{i1}(in_gate~=1,:);
            in_gate01 = SigmaGating_com(XY_1,mu(1,:),Sigma(:,:,1),model);
            in_gate02 = SigmaGating_com(XY_2,mu(2,:),Sigma(:,:,2),model);
            %α�Ӳ�
            C_f = [C_f;XY_1(in_gate01~=1,:);XY_2(in_gate02~=1,:)];
            %�ɴ��ݵ�����
            XY_1 = XY_1(in_gate01==1,:);
            XY_2 = XY_2(in_gate02==1,:);

            cell_spliting{idi_emp(1)} = XY_1;
            idi_emp(1) = [];
            cell_spliting{idi_emp(1)} = XY_2;
            idi_emp(1) = [];
            cell_spliting{i1} = [];
        end
    end
end
idi0 = zeros(N_splitMax,1);
for i0 = 1:N_splitMax %�ҳ������鼯������
   if ~isempty(cell_spliting{i0})
        idi0(i0) = 1;
    end
end
if any(idi0)
    idi = find(idi0==1);
    for i = idi
        C_f = [C_f;cell_spliting{i}];
    end  
end
%�Ӳ����·���
[cell_meas,C_t] = ClutterRedistribute_com(cell_splited,N_splited-1,C_f,model);

N_cell_meas = size(cell_meas,1);
ksem_Z_idi = cell(N_cell_meas,1);
for i = 1:N_cell_meas
    idi = find(sum(ismember(W,cell_meas{i}),2))';
    ksem_Z_idi{i} = idi;
end
ksem_Z = cell_meas;
end