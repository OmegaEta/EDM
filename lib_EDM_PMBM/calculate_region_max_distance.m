function [max_distances, region_angles] = calculate_region_max_distance(C, z_bar, num_regions)
    % 输入参数：
    % C: 2×N矩阵，表示所有点的笛卡尔坐标（每列为一个点）
    % z_bar: 2×1向量，质心坐标
    % num_regions: 分区域数量（例如36表示每10度一个区域）
    % 输出参数：
    % max_distances: 每个区域的最远投影距离
    % region_angles: 每个区域的中心角度（弧度）

    % 步骤1：计算所有点相对于质心的极坐标
    temp = C - z_bar; % 2×N，笛卡尔坐标差
    [theta, r] = cart2pol(temp(1,:), temp(2,:)); % theta∈[-π, π]
    theta = mod(theta, 2*pi); % 转换到[0, 2π)

    % 步骤2：划分角度区域
    edges = linspace(0, 2*pi, num_regions + 1);
    region_angles = (edges(1:end-1) + edges(2:end)) / 2; % 区域中心角度
    [~, bin] = histc(theta, edges);
    bin(bin == num_regions + 1) = 1; % 处理边界点

    % 步骤3：预计算所有点的垂直方向投影向量
    perp_angles = region_angles(bin) + pi/2; % 每个点所属区域的垂直方向
    % 笛卡尔单位向量计算
    unit_vectors = [cos(perp_angles); sin(perp_angles)]; % 2×N

    % 步骤4：计算投影距离（绝对值）
    projection_dist = sum(temp .* unit_vectors, 1); % 点积计算投影
    projection_dist = abs(projection_dist); % 取绝对值表示距离

    % 步骤5：按区域求最大值
    max_distances = accumarray(bin(:), projection_dist(:), ...
                               [num_regions, 1], @max, 0);

    % 处理空区域（可选：用NaN标记）
    region_counts = accumarray(bin(:), 1);
    max_distances(region_counts == 0) = NaN;
end