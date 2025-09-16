function in_gate = irregularGating(W,GGIW,model)

%irregular
if isempty(W)
    in_gate = [];
    return;
end

d = 2;
in_gate = false(size(W,2),1);
H = model.measmodel.H(GGIW.m);

% if all(GGIW.shape.m(3:4)==0)
if 1
    %原方法
    S = 0.5*GGIW.ScaleEplision/(GGIW.v-2*d-2) + H*GGIW.P*H' + model.measmodel.R;
    S = (S + S')/2;
    
    nu = W - repmat(model.measmodel.h(GGIW.m),[1,size(W,2)]);
    dist= sum((inv(chol(S))'*nu).^2);
    
    in_gate(dist<model.gamma) = true;
else
    S = GGIW.shape.shape_parameter./reshape(GGIW.shape.degrees_freedom-2*d-2,1,1,model.direction) + H*GGIW.P*H' + model.measmodel.R;
    S = (S + permute(S, [2 1 3]))/2;
    
    %nu为每个量测到质心的距离
    nu = W - repmat(model.measmodel.h(GGIW.shape.m),[1,size(W,2)]);
    
%     measureinblock_index = irregular_PerBlockCov_gating(nu,model);

    % 预分配存储逆Cholesky因子的三维数组
    invL = zeros(2, 2, model.direction);
    dist = zeros(1, model.direction);
    for i = 1:model.direction
        L = chol(S(:,:,i), 'lower');   % Cholesky分解得到下三角矩阵L
        invL(:,:,i) = inv(L);          % 计算L的逆矩阵

        transformed = invL(:,:,i)' * nu(:,i);  % 每个向量独立计算
        dist(i) = sum(transformed.^2);
    end

    in_gate(dist < model.gamma) = true;
end



end
