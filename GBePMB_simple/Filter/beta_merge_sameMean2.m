function [alpha_k,beta_k] = beta_merge_sameMean2(alphas,betas,ws)
%
N = length(ws);

mu_w = (alphas./(alphas+betas))*ws;
beta_k = betas*ws;
beta_k = 1;
alpha_k = (mu_w/(1-mu_w))*beta_k;


% alpha_k = beta_k;
% alpha_k = 10;
% beta_k = 10;
C1 = psi(alpha_k);
C2 = psi(beta_k);
C3 = psi(alpha_k+beta_k);

A = zeros(N,2);
% A1 = zeros(N,2);
B = zeros(N,1);
% x=0.00001:0.0001:0.99999;
for i = 1:N
    % y=(x.^(alphas(i)-1)).*((1-x).^(betas(i)-1)).*log(x);
    % A(i,1) = trapz(x,y);
    % 
    % y=(x.^(alphas(i)-1)).*((1-x).^(betas(i)-1)).*log(1-x);
    % A(i,2) = trapz(x,y);

    fun1 = @(x) (x.^(alphas(i)-1)).*((1-x).^(betas(i)-1)).*log(x);
    fun2 = @(x) (x.^(alphas(i)-1)).*((1-x).^(betas(i)-1)).*log(1-x);
    A(i,1) = integral(fun1,0,1);
    A(i,2) = integral(fun2,0,1);

    % B(i) = ws(i)*gamma(alphas(i)+betas(i))/(gamma(alphas(i))*gamma(betas(i)));
    B(i) = ws(i)*(1/beta(alphas(i),betas(i)));
end
D = A'*B;
iter = 0;

i_s = 0;
S = [1 0.5 0.1 0.05 0.01 0.002];
i_s = i_s+1;
step = S(i_s);

CE_k = (alpha_k-1)*D(1) + (beta_k-1)*D(2) - log(gamma(alpha_k)) - log(gamma(beta_k)) + log(gamma(alpha_k+beta_k));

fprintf("\n-------beta_merge_sameMean2------\n");
while iter<1000
    iter = iter + 1;

    beta_steer = (mu_w/(1-mu_w))*D(1) + D(2) - (mu_w/(1-mu_w))*C1 - C2 + (1/(1-mu_w))*C3;
    
    % alpha_new = alpha_k - alpha_steer/ab;
    % beta_new = beta_k - beta_steer/ab;

    beta_new = beta_k +beta_steer*step;
    alpha_new = (mu_w/(1-mu_w))*beta_new;
    
    CE_new = (alpha_new-1)*D(1) + (beta_new-1)*D(2) - log(gamma(alpha_new)) - log(gamma(beta_new)) + log(gamma(alpha_new+beta_new));

    CE_c = (CE_new-CE_k)/CE_k;

    
    % abs(KLDc)<0.0002 || 
    % (abs(beta_steer)<0.97 &&
    if  abs(CE_c)<0.00005
        fprintf("\n%d:CE_k=%f;beta=%f;\nbeta_steer=%f;CE_c=%f \n",iter,CE_k,beta_k,beta_steer,CE_c);
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

        pd = makedist('Beta','a',alpha_k,'b',beta_k);
        h = pd.plot;

        % fprintf("\n%d:CE_k=%f;beta=%f;\nbeta_steer=%f;CE_c=%f \n",iter,CE_k,beta_k,beta_steer,CE_c);
        delete(h);
        % alpha_k = alpha_new;
        % beta_k = beta_new;
        
        C1 = psi(alpha_k);
        C2 = psi(beta_k);
        C3 = psi(alpha_k+beta_k);
    end

end
fprintf("\nbeta_merge_sameMean2\n%d:CE_k=%f;beta=%f;\nbeta_steer=%f;CE_c=%f \n\n",iter,CE_k,beta_k,beta_steer,CE_c);
fprintf("\n========beta_merge_sameMean2========\n");
end