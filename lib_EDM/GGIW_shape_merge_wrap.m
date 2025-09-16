function [w_hat,GGIW_hat] = GGIW_shape_merge_wrap(w,GGIWs)

direction = length(GGIWs(end).num_parameter);
shape_array = zeros(2,2,direction,length(w));
for i = 1:length(w)
    shape_array(:,:,:,i) = cat(3, GGIWs(i).shape_parameter);
end
 
[w_hat,GGIW_hat.num_parameter,GGIW_hat.inverse_scale,GGIW_hat.m,GGIW_hat.deviation_matrix,GGIW_hat.degrees_freedom,GGIW_hat.shape_parameter] = ...
    GGIW_shapemerge(exp(w), reshape([GGIWs.num_parameter],[1,direction,length(w)]),...
     reshape([GGIWs.inverse_scale],[1,direction,length(w)]),[GGIWs.m],...
    reshape([GGIWs.deviation_matrix],[4,4,length(w)]),reshape([GGIWs.degrees_freedom],[1,direction,length(w)])...
    ,shape_array);

end

