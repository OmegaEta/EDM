function in_gate = SigmaGating_com(W,m,Sigma,model)
% clarify
if isempty(W)
    in_gate = [];
    return;
end

d = 2;

in_gate = false(size(W,1),1);

S = Sigma;
nu = W - repmat(m,[size(W,1),1]);
dist= sum((inv(chol(S))'*nu').^2);

in_gate(dist<model.gamma) = true;

end
