function nc = optimizeClutter(ng_in,psg,lambda_X,Tc)
% 最优化计算杂波估计
if sum(psg)<=0 || sum(ng_in)<=0
    nc = 0;
else
    if Tc == 1
        nc = (ng_in-lambda_X)/psg;
        % if nc<0
        %     nc = 0;
        % end
    
    else
        fun = @(x)(sum((ng_in.*psg)./(psg.*x+lambda_X))-sum(psg))^2;
    
        [nc,fval] = fminunc(fun,60);
    end

end

end

