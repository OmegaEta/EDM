function area = GetArea(GGIW,model)
%
d = 2;
if GGIW.v==0
    area = 0;
    return 
end
H = model.measmodel.H(GGIW(end).m);
Sigma2_ = GGIW.V/(GGIW.v-2*d-2) + model.measmodel.R;
Sigma2_ = 0.5*(Sigma2_+Sigma2_');
mu_ = [GGIW.m(1) GGIW.m(2)];

[~,Sigma] = eig(Sigma2_);
a = sqrt(Sigma(1));
b = sqrt(Sigma(end));

area = pi*a*b;
end

