function [alpha_k,beta_k] = beta_merge(alphas,betas,ws)
%
N = length(ws);
alpha_k = alphas*ws;
beta_k = betas*ws;

% alpha_k = beta_k;
% alpha_k = 10;
% beta_k = 10;
C1 = psi(alpha_k);
C2 = psi(beta_k);
C3 = psi(alpha_k+beta_k);

A = zeros(N,2);
B = zeros(N,1);
x=0.00001:0.01:0.99999;
for i = 1:N
    y=(x.^(alphas(i)-1)).*((1-x).^(betas(i)-1)).*log(x);
    A(i,1) = trapz(x,y);

    y=(x.^(alphas(i)-1)).*((1-x).^(betas(i)-1)).*log(1-x);
    A(i,2) = trapz(x,y);
    
    B(i) = ws(i)*gamma(alphas(i)+betas(i))/(gamma(alphas(i))*gamma(betas(i)));
end
D = A'*B;
iter = 1;
step = 0.2;%0.05
while iter<10000
    iter = iter + 1;
    alpha_steer = D(1) - C1 + C3;
    beta_steer = D(2) - C2 + C3;
    ab = (alpha_k-1)*D(1) + (beta_k-1)*D(2) - log(gamma(alpha_k)) - log(gamma(beta_k)) + log(gamma(alpha_k+beta_k));
    % alpha_new = alpha_k - alpha_steer/ab;
    % beta_new = beta_k - beta_steer/ab;
    alpha_new = alpha_k +alpha_steer*step;
    beta_new = beta_k +beta_steer*step;
    
    err = alpha_new-alpha_k;
    % if mod(iter,1000)==1
    %     fprintf("%d:%f\n",iter,abs(err));
    % end
    
    if abs(err)<0.00001
        break;
    else
        if beta_new<1 && alpha_new<1
            break;
        elseif alpha_new<1
            beta_k = beta_new;
        elseif beta_new<1
            alpha_k = alpha_new;
        else
            alpha_k = alpha_new;
            beta_k = beta_new;
        end
        % alpha_k = alpha_new;
        % beta_k = beta_new;
        C1 = psi(alpha_k);
        C2 = psi(beta_k);
        C3 = psi(alpha_k+beta_k);
    end

end
fprintf("%d:%f;alpha=%f;beta=%f\n",iter,abs(err),alpha_k,beta_k)            ;
end