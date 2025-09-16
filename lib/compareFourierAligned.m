function [min_diff, phi_star] = compareFourierAligned(A, B)
%     % 采样角度
%     theta = linspace(0, 2*pi, 360);  % 更高精度角度分辨率
%     dtheta = theta(2) - theta(1);    % 角度间隔
%     
%     % 枚举相位平移 φ
%     phi_vals = linspace(0, 2*pi, 360);  % 相位扫描
%     diffs = zeros(size(phi_vals));
%     
%     for i = 1:length(phi_vals)
%         phi = phi_vals(i);
%         
%         % 将A进行相位平移：evaluate at (theta + phi) mod 2π
%         shifted_theta = mod(theta + phi, 2*pi);
%         A_shifted = evaluateFourierSeries(A.a, A.b, shifted_theta);
%         
%         % 原始B
%         B_vals = evaluateFourierSeries(B.a, B.b, theta);
%         
%         % 计算L2差异
%         diffs(i) = sqrt(sum((A_shifted - B_vals).^2) * dtheta);
%     end
%     
%     % 找到最小差异对应的phi
%     [min_diff, idx] = min(diffs);
%     phi_star = phi_vals(idx);
%     
% %     % 可选：绘图查看匹配效果
% %     figure;
% %     plot(theta, evaluateFourierSeries(A.a, A.b, mod(theta + phi_star, 2*pi)), 'r', 'LineWidth', 2); hold on;
% %     plot(theta, evaluateFourierSeries(B.a, B.b, theta), 'b--', 'LineWidth', 2);
% %     legend('Aligned A', 'B');
% %     title(['Aligned with \phi = ', num2str(phi_star)]);

%%
% 高密度角度采样
    theta = linspace(0, 2*pi, 360).';   % [360×1]
    dtheta = theta(2) - theta(1);

%     % 枚举所有平移角 phi
%     phi_vals = linspace(0, 2*pi, 360); % [1×360]
%     
%     % 构建 θ + φ 的矩阵，大小为 [360 × 360]
%     Theta_mat = mod(theta + phi_vals, 2*pi);  % [360×360]

    % 批量评估 A(θ + φ)
    A_shifted_all = evaluateFourierSeries(2*A.a, 2*A.b, theta);  % [360×360]

    % 评估 B(θ)，只需一次
    B_vals = evaluateFourierSeries(B.a, B.b, theta);  % [360×1]

    % 广播相减，计算每一列误差
    diff_mat = A_shifted_all - B_vals;  % 自动广播 [360×360]

    % 每列计算 L2 范数（积分近似）
    diffs = sqrt(sum(diff_mat.^2, 1) * dtheta);  % [1×360]

    % 最小值和对应φ
    [min_diff, idx] = min(diffs);
%     phi_star = phi_vals(idx);
    phi_star = 0;

end
