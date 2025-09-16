function Z = productRandomMeas(TrackLocation,lamda)
% 根据运动轨迹位置 产生量测

N = size(TrackLocation,1);
Z = cell(N,1);

for i = 1:N
    count = random('Poisson',lamda);
    xX = TrackLocation{i};
    x = xX.x;
    X = xX.X;
    z = mvnrnd(x(1:2),X,count);
    Z{i} = z;
    
end
end

