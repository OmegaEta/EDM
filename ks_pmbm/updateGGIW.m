function [GGIW,lik] = updateGGIW(GGIW,W,model)

d = 2;

card_W = size(W,2);
if card_W==0 
%     GGIW = struct('a',0,'b',1,'m',[0;0],'v',0,'P',zeros(4,4),'V',zeros(2,2));
    lik = 0;
    return;
end

GGIW_.a = GGIW(end).a + card_W;
GGIW_.b = GGIW(end).b + 1;

z_bar = mean(W,2);
epsilon = z_bar - model.measmodel.h(GGIW(end).m);
H = model.measmodel.H(GGIW(end).m);

X_hat = GGIW(end).V/(GGIW(end).v - 2*d - 2);
X_hat = (X_hat + X_hat')/2;

S = H*GGIW(end).P*H' + X_hat/card_W;
S = (S + S')/2;

Vs= chol(S); det_S= prod(diag(Vs))^2; inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';

K = GGIW(end).P*H'*iS;

GGIW_.m = GGIW(end).m + K*epsilon;
GGIW_.P = GGIW(end).P - K*H*GGIW(end).P;

temp = (W - z_bar);
Z = temp*temp';

X_sqrt = sqrtm_2by2(X_hat); %sqrtm?
S_sqrt_inv = sqrtm(iS);
N = X_sqrt*S_sqrt_inv*(epsilon*epsilon')*S_sqrt_inv'*X_sqrt';

GGIW_.v = GGIW(end).v + card_W;
GGIW_.V = GGIW(end).V + N + Z;

%

lik = (GGIW(end).v-d-1)/2*log(det2(GGIW(end).V)) - (GGIW_.v-d-1)/2*log(det2(GGIW_.V))...
        + gamma2ln((GGIW_.v-d-1)/2) - gamma2ln((GGIW(end).v-d-1)/2)...
        + log(det2(X_hat))/2 - log(det_S)/2 + gammaln(GGIW_.a)...
        - gammaln(GGIW(end).a) + GGIW(end).a*log(GGIW(end).b)  - GGIW_.a*log(GGIW_.b)...
        - (card_W*log(pi)+log(card_W))*d/2;
    
GGIW(end) = GGIW_;

% X_hat = GGIW(end).V/(GGIW(end).v - 2*d - 2);
% X_hat = (X_hat + X_hat')/2;
% Sigmacircle(GGIW(end).m(1),GGIW(end).m(2),X_hat,2,7);

end

