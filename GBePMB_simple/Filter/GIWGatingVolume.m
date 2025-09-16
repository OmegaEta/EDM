function V = GIWGatingVolume(Cr,v,V,paras,measmodel)
% 计算GIW门限内的面积

S = measmodel.H*Cr*measmodel.H' + ...
    measmodel.Ch*V/(v-6) + measmodel.Cv;
S = (S+S')/2;


V = ellipsoidalVolume(S,paras.gating.size,measmodel);
end

