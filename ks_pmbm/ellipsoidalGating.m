function in_gate = ellipsoidalGating(W,GGIW,model)

if isempty(W)
    in_gate = [];
    return;
end

d = 2;

in_gate = false(size(W,2),1);

H = model.measmodel.H(GGIW.m);
S = GGIW.V/(GGIW.v-2*d-2) + H*GGIW.P*H' + model.measmodel.R;
S = (S + S')/2;
nu = W - repmat(model.measmodel.h(GGIW.m),[1,size(W,2)]);
dist= sum((inv(chol(S))'*nu).^2);

in_gate(dist<model.gamma) = true;
% 
S0 = sqrtm(S);

% figure(1);
% Sigmacircle(GGIW.m(1),GGIW.m(2),S,2,6);
%legend('S gating');
% W = W(:,in_gate);
% plot(W(1,:),W(2,:),'go');
% 
% S0=1.5*3.5*S;
% Sigmacircle(GGIW.m(1),GGIW.m(2),S0,2,5);
% %legend('S*2.25 gating');

end
