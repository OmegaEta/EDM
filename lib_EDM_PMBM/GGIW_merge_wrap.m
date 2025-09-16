function [w_hat,GGIW_hat] = GGIW_merge_wrap(w,GGIWs)
% % % %测试用例
% w = [w;-14.47];
% GGIWs(2).a = GGIWs.a;
% GGIWs(2).b = GGIWs.b;
% GGIWs(2).m = GGIWs(1).m.*1.05;
% GGIWs(2).P = GGIWs(1).P;
% GGIWs(2).v = GGIWs.v;
% GGIWs(2).V = GGIWs.V;
% 
% GGIWs(1).V = [1,2;3,4];
% GGIWs(2).V = [11,21;31,41];
% ans(1,:,2)
% for i = 1:length(w)
%     GGIWs.V
% end

[w_hat,GGIW_hat.a,GGIW_hat.b,GGIW_hat.m,GGIW_hat.P,GGIW_hat.v,GGIW_hat.V,GGIW_hat.ScaleEplision] = ...
    GGIW_merge(exp(w),[GGIWs.a],[GGIWs.b],[GGIWs.m],...
    reshape([GGIWs.P],[4,4,length(w)]),[GGIWs.v],reshape([GGIWs.V],[2,2,length(w)]),...
    reshape([GGIWs.ScaleEplision],[2,2,length(w)]));

GGIW_hat.shape = GGIWs.shape;
GGIW_hat.Shape_coefficients = GGIWs.Shape_coefficients;

end

