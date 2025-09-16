function adaGauss = adaMeas(XY,model)
%clarify
N_XY = size(XY,1);
if N_XY==0
    XY_bar = [1000 1000];
    Z = [0 0;0 0];
else
    XY_bar = mean(XY,1);
    temp = (XY - XY_bar);
    Z = temp'*temp;
end
adaGauss.m = XY_bar;
adaGauss.P = Z+model.ada.noise;
% Z1 = sum(temp.*temp,2);
% GGIW0.a = 700;
% GGIW0.b = 100;
% GGIW0.m = [XY_bar 0 0]';
% GGIW0.P = diag([1 1 0 0]);
% GGIW0.V = Z + [1 0;0 1];%model.deri.P;
% GGIW0.v = N_XY + 15;
% 
% [GGIW,lik] = updateGGIW(GGIW0,XY',model);
end

