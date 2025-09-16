function [estimates,trajectoryEstimates] = estimator_irregular(MBM,model)

trajectoryEstimates = [];
num_edge = model.direction*10;
estimates.g = cell(1,0);%gamma
estimates.x = zeros(4,0);%Shape component Difference.
estimates.X = zeros(2,2,0);%Shape component Extension covariance.
estimates.r = cell(0,1);%Target shape parameters Sc.
estimates.shape = zeros(2,2,model.direction,0);%Shape covariance of components.

d = 2;
% edge_angle = linspace(0,2*pi,num_edge);
[~,idx] = max(MBM.w);
table_entry = MBM.table(idx,:);
for i = 1:length(table_entry)
    if table_entry(i) > 0 && MBM.track{i}(table_entry(i)).Bern.r >= model.exist_r
        [~,ideath] = max(MBM.track{i}(table_entry(i)).Bern.w_death);
        if ideath == length(MBM.track{i}(table_entry(i)).Bern.w_death)
            g = MBM.track{i}(table_entry(i)).Bern.GGIW(end).shape.num_parameter./MBM.track{i}(table_entry(i)).Bern.GGIW(end).shape.inverse_scale;
            estimates.g = cat(2,estimates.g,g);
            estimates.x = [estimates.x MBM.track{i}(table_entry(i)).Bern.GGIW(end).shape.m];

            r.a = MBM.track{i}(table_entry(i)).Bern.GGIW(end).Shape_coefficients.a;
            r.b = MBM.track{i}(table_entry(i)).Bern.GGIW(end).Shape_coefficients.b;
            r_ = evaluateFourierSeries(r.a,r.b,model.direction_angle);
            [x,y] = pol2cart(model.direction_angle,r_);range_ = [x;y];
            temp = range_-mean(range_,2);
            X_ = (temp*temp'/size(temp,2));
            X = X_;
            estimates.X = cat(3,estimates.X,X);
            estimates.r = cat(1,estimates.r,r);

            r_ = evaluateFourierSeries(r.a,r.b,model.direction_angle);
            [x, y] = pol2cart(model.direction_angle, r_);
            range_ = [x;y];
            A = range_-mean(range_,2);
            Cov_irregular(:,:,:) = (reshape(A, 2, 1, []) .* permute(reshape(A, 2, 1, []), [2, 1, 3]));
            shape = Cov_irregular;
            estimates.shape = cat(4,estimates.shape,shape);
        end
        t_death = MBM.track{i}(table_entry(i)).Bern.t_death(ideath);
        tlen = t_death - MBM.track{i}(table_entry(i)).Bern.t_birth + 1;

        trajectoryEstimates(end+1,1).t_birth = [MBM.track{i}(table_entry(i)).assocHistory(1).t];
        trajectoryEstimates(end,1).t_death = t_death;
        trajectoryEstimates(end,1).g = [MBM.track{i}(table_entry(i)).Bern.GGIW(1:tlen).a]./[MBM.track{i}(table_entry(i)).Bern.GGIW(1:tlen).b];
        trajectoryEstimates(end,1).x = [MBM.track{i}(table_entry(i)).Bern.GGIW(1:tlen).m];
        
        temp = arrayfun(@(x) x.ScaleEplision/(x.v-2*d-2), MBM.track{i}(table_entry(i)).Bern.GGIW(1:tlen), 'un', false);
        trajectoryEstimates(end,1).X = zeros(d,d,length(temp));
        for j = 1:length(temp)
            trajectoryEstimates(end,1).X(:,:,j) = temp{j};
        end
        
    end
%     trajectoryEstimates.X = cat(4,trajectoryEstimates.X,X);
end  
    
    
end

