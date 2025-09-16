%更新膨胀函数，赋予GGIW膨胀参数
function tracks = ExtensionDifferenceFunc(track,W,model)

% boundaryMapp = cell(length(track));
%计算两个形状的相交
for i = 1:length(track)
    for j = 1:length(track{i,1})
        internal_points = W(:,track{i,1}(j).assocHistory(end).meas);
    %     original_boundary_ =  evaluate_Fourier(track{i,1}.Bern.GGIW(end).Shape_coefficients, model.direction_angle);
        
        original_boundary_ = evaluateFourierSeries(track{i,1}(j).Bern.GGIW(end).Shape_coefficients.a,track{i,1}(j).Bern.GGIW(end).Shape_coefficients.b,model.direction_angle);
        x = original_boundary_ .* cos(model.direction_angle)+ track{i}(j).Bern.GGIW(end).m(1); % X坐标
        y = original_boundary_ .* sin(model.direction_angle)+ track{i}(j).Bern.GGIW(end).m(2); % Y坐标
        original_boundary = [x;y];
        
        %膨胀函数
        %尺度椭圆
        [Eplision.V,Eplision.D] = eig(2*track{i,1}(j).Bern.GGIW(end).ScaleEplision/(track{i,1}(j).Bern.GGIW(end).v-2d-2));
        Eplision.center = track{i,1}(j).Bern.GGIW(end).m(1:2)';
%         Eplision.center = mean(original_boundary,2)';

        % 执行膨胀变形
        [coefficients_EDM] = shape_expansion(...
        original_boundary', Eplision.V, Eplision.D, Eplision.center);
    
        %更新膨胀函数
        track{i, 1}(j).Bern.GGIW(end).shape.EDM_coefficients = coefficients_EDM;

        %
        tracks = track;
    end
end
end