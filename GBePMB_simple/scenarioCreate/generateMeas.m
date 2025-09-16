function Z = generateMeas(gt,measmodel,motionmodel,T)
% 产生量测


nt = size(gt,2);
Z = cell(T,1);
%target-generated measurements
for i = 1:nt
    for t = gt(i).x_bt:gt(i).x_dt
        %generate measurements only if the target is detected
        %no misdetection at time 1
        if rand < measmodel.Pd || t == gt(i).x_bt
            nz = poissrnd(gt(i).x_lambda(t-gt(i).x_bt+1));
            Z{t} = [Z{t} mvnrnd(gt(i).xr(1:2,t-gt(i).x_bt+1)',...
                measmodel.Ch*gt(i).X(:,:,t-gt(i).x_bt+1)+measmodel.Cv,nz)'];
        end
    end
end

%append clutter
for i = 1:T
    nz = poissrnd(measmodel.c_lambda);
    zx = rand(nz,1)*(motionmodel.range_x(2)-motionmodel.range_x(1)) + motionmodel.range_x(1);
    zy = rand(nz,1)*(motionmodel.range_y(2)-motionmodel.range_y(1)) + motionmodel.range_y(1);
    Z{i} = [Z{i} [zx zy]'];
end
end

