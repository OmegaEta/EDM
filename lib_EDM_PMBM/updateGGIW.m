function [GGIW,lik] = updateGGIW(GGIW,W,model)

%Apply EDM
new_internal = internal_expansion(W,GGIW(end).shape.EDM_coefficients);
W_ = W;
%Change to the new measurement set after mapping.
W = new_internal';

d = 2;
card_W = size(W,2);

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

X_sqrt = sqrtm_2by2(X_hat);
S_sqrt_inv = sqrtm_2by2(iS);
% S_sqrt_inv = sqrtm(iS);
N = X_sqrt*S_sqrt_inv*(epsilon*epsilon')*S_sqrt_inv'*X_sqrt';

GGIW_.v = GGIW(end).v + card_W;
GGIW_.V = GGIW(end).V + N + Z;

%形状估计
% GGIW_.shape = GGIW.shape;

lik = (GGIW(end).v-d-1)/2*log(det2(GGIW(end).V)) - (GGIW_.v-d-1)/2*log(det2(GGIW_.V))...
        + gamma2ln((GGIW_.v-d-1)/2) - gamma2ln((GGIW(end).v-d-1)/2)...
        + log(det2(X_hat))/2 - log(det_S)/2 + gammaln(GGIW_.a)...
        - gammaln(GGIW(end).a) + GGIW(end).a*log(GGIW(end).b) - GGIW_.a*log(GGIW_.b)...
        - (card_W*log(pi)+log(card_W))*d/2;

% GGIW(end).shape=predictShapeRotate(GGIW(end).m(3:4),GGIW_.m(3:4),GGIW(end).shape);

%传输参数
GGIW_.shape = GGIW(end).shape;
GGIW_.Shape_coefficients = GGIW(end).Shape_coefficients;

%尺度椭圆
%尺度椭圆的更新
epsilon = z_bar - model.measmodel.h(GGIW(end).m);

X_hat_ = GGIW(end).ScaleEplision/(GGIW(end).v - 2*d - 2);
X_hat_ = (X_hat_ + X_hat_')/2;

S_ = H*GGIW(end).P*H' + X_hat_/card_W;
S_ = (S_ + S_')/2;

Vs= chol(S_);inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';

temp = (W_ - z_bar);
Z = temp*temp';

X_sqrt = sqrtm_2by2(X_hat);
S_sqrt_inv = sqrtm_2by2(iS);
% S_sqrt_inv = sqrtm(iS);
N = X_sqrt*S_sqrt_inv*(epsilon*epsilon')*S_sqrt_inv'*X_sqrt';

GGIW_.ScaleEplision = GGIW(end).ScaleEplision + N + Z;

GGIW(end) = GGIW_;

end

