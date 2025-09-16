function [alpha_k,beta_k] = beta_merge_nobias(states,w)
%
w = exp(w);

idx = w > 0;
s = [states.s];
t = [states.t];
ws = w(idx);
alphas = s(idx);
betas = t(idx);


N = length(ws);
if N<=1
    alpha_k = alphas;
    beta_k = betas;
    return;
end

mu_w = (alphas./(alphas+betas))*ws;
beta_k = betas*ws;
alpha_k = (mu_w/(1-mu_w))*beta_k;
if beta_k<1
   beta_k = 1.01;
   alpha_k = (mu_w/(1-mu_w))*beta_k;
end
C1 = psi(alpha_k);
C2 = psi(beta_k);
C3 = psi(alpha_k+beta_k);

A = zeros(N,2);
B = zeros(N,1);
for i = 1:N
    fun1 = @(x) (x.^(alphas(i)-1)).*((1-x).^(betas(i)-1)).*log(x);
    fun2 = @(x) (x.^(alphas(i)-1)).*((1-x).^(betas(i)-1)).*log(1-x);
    A(i,1) = integral(fun1,0,1);
    A(i,2) = integral(fun2,0,1);

    B(i) = ws(i)*(1/beta(alphas(i),betas(i)));
end
D = A'*B;
iter = 0;

i_s = 0;
S = [1 0.5 0.1 0.05 0.01 0.002];
i_s = i_s+1;
step = S(i_s);

CE_k = (alpha_k-1)*D(1) + (beta_k-1)*D(2) - log(gamma(alpha_k)) - log(gamma(beta_k)) + log(gamma(alpha_k+beta_k));

while iter<1000
    iter = iter + 1;

    beta_steer = (mu_w/(1-mu_w))*D(1) + D(2) - (mu_w/(1-mu_w))*C1 - C2 + (1/(1-mu_w))*C3;
    
    beta_new = beta_k +beta_steer*step;
    alpha_new = (mu_w/(1-mu_w))*beta_new;
    
    CE_new = (alpha_new-1)*D(1) + (beta_new-1)*D(2) - log(gamma(alpha_new)) - log(gamma(beta_new)) + log(gamma(alpha_new+beta_new));
    CE_c = (CE_new-CE_k)/CE_k;

    if  abs(CE_c)<0.00005
        break;
    else
        if (beta_new<=1 || alpha_new<=1)
            if i_s < length(S)
                i_s = i_s + 1;
                step = S(i_s);
                continue;
            else
               break;
            end
        end
        alpha_k = alpha_new;
        beta_k = beta_new;
        CE_k = CE_new;
    
        C1 = psi(alpha_k);
        C2 = psi(beta_k);
        C3 = psi(alpha_k+beta_k);
    end

end
end