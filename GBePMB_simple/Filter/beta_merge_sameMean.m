function [alpha_k,beta_k] = beta_merge_sameMean(alphas,betas,ws)
%
N = length(ws);

weight_mean = (alphas./(alphas+betas))*ws;
beta_k = betas*ws;
% beta_k = 5;
alpha_k = (weight_mean/(1-weight_mean))*beta_k;


% alpha_k = beta_k;
% alpha_k = 10;
% beta_k = 10;
C1 = psi(alpha_k);
C2 = psi(beta_k);
C3 = psi(alpha_k+beta_k);

A = zeros(N,2);
A1 = zeros(N,2);
B = zeros(N,1);
x=0.00001:0.0001:0.99999;
for i = 1:N
    % y=(x.^(alphas(i)-1)).*((1-x).^(betas(i)-1)).*log(x);
    % A(i,1) = trapz(x,y);
    % 
    % y=(x.^(alphas(i)-1)).*((1-x).^(betas(i)-1)).*log(1-x);
    % A(i,2) = trapz(x,y);

    fun1 = @(x) (x.^(alphas(i)-1)).*((1-x).^(betas(i)-1)).*log(x);
    fun2 = @(x) (x.^(alphas(i)-1)).*((1-x).^(betas(i)-1)).*log(1-x);
    A1(i,1) = integral(fun1,0,1);
    A1(i,2) = integral(fun2,0,1);

    % B(i) = ws(i)*gamma(alphas(i)+betas(i))/(gamma(alphas(i))*gamma(betas(i)));
    B(i) = ws(i)*(1/beta(alphas(i),betas(i)));
end
D = A1'*B;
iter = 0;
step = 0.01;%0.05
CE_k = (alpha_k-1)*D(1) + (beta_k-1)*D(2) - log(gamma(alpha_k)) - log(gamma(beta_k)) + log(gamma(alpha_k+beta_k));

fprintf("\n---------beta_merge_sameMean---------\n");

while iter<1000
    iter = iter + 1;

    beta_steer = (weight_mean/(1-weight_mean))*D(1) + D(2) - C2 + C3;
    
    % alpha_new = alpha_k - alpha_steer/ab;
    % beta_new = beta_k - beta_steer/ab;

    beta_new = beta_k +beta_steer*step;
    alpha_new = (weight_mean/(1-weight_mean))*beta_new;
    % CE_new = (alpha_k-1)*D(1) + (beta_k-1)*D(2) - log(gamma(alpha_k)) - log(gamma(beta_k)) + log(gamma(alpha_k+beta_k));
    CE_new = (alpha_new-1)*D(1) + (beta_new-1)*D(2) - log(gamma(alpha_new)) - log(gamma(beta_new)) + log(gamma(alpha_new+beta_new));

    % err = (beta_new-beta_k)/beta_k;

    CE_c = (CE_new-CE_k)/CE_k;
    % if mod(iter,1000)==1
    %     fprintf("%d:%f\n",iter,abs(err));
    % end

    
    % abs(KLDc)<0.0002 || 
    % (abs(beta_steer)<0.97 &&
    if  abs(CE_c)<0.00005
        fprintf("\n%d:CE_k=%f;beta=%f;\nbeta_steer=%f;CE_c=%f \n",iter,CE_k,beta_k,beta_steer,CE_c);
        break;
    else
        if beta_new<=1 || alpha_new<=1
            break;
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
fprintf("\nbeta_merge_sameMean\n%d:CE_k=%f;beta=%f;\nbeta_steer=%f;CE_c=%f \n\n",iter,CE_k,beta_k,beta_steer,CE_c);
fprintf("\n========beta_merge_sameMean1========\n");
end