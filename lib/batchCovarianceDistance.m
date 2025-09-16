function D = batchCovarianceDistance(P1, P2)
    % P1, P2: 2×2×36×n
    % D(i,j): shape distance for i-th component of j-th target
    [~, ~, K, N] = size(P1);
    D = zeros(K, N);

    for j = 1:N
        for i = 1:K
            A = P1(:,:,i,j);
            B = P2(:,:,i,j);
            sqrtA = sqrtm_2by2(A);  % 替代 sqrtm(A)，速度更快
            inner = sqrtA * B * sqrtA;
            sqrtInner = sqrtm_2by2(inner);
            D(i,j) = trace(A + B - 2 * sqrtInner)/2;
        end
    end
end
