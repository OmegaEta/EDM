% The structure of V is a 2x2 matrix for each direction and each GGIW, so the structure is 2x2xdirectionxnum(GGIW).
function [w2,a2,b2,m2,P2,nu2,V2] = GGIW_shapemerge(w,a,b,m,P,nu,V)

nu_min = 7;
v_min = 1;

% Dimensions
n_x = size(m,1);
d = size(V,1);
% Number of components
N=length(w);

% Sum of weights
wb = sum(w);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct merging components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c1 = 0;
c2 = 0;
C1 = zeros(d,d);
C2 = 0;
C3 = 0;
nu_ = reshape(nu, 1, 1, length(a), []);

for i = 1:N
    % Gammas
    c1 = c1+w(i)*(psi(0,a(:,:,i))-log(b(:,:,i)));
    c2 = c2+w(i)*a(:,:,i)./b(:,:,i);
    % inverse Wisharts  
    C1 = C1+w(i)*(nu_(:,:,:,i)-d-1).*(V(:,:,:,i).\eye(d));
    C2 = C2+w(i)*sum(psi(0,(nu(:,:,i)-(d+(1:d)'))./2));
    C3 = C3 + w(i)*log((V(1,1,:,i).*V(2,2,:,i) - V(1,2,:,i).*V(2,1,:,i)));
end
C3 = reshape(C3,1,[]);
c = c1./wb - log(c2./wb);
C = d*wb*log(wb) -wb*reshape(log((C1(1,1,:).*C1(2,2,:) - C1(1,2,:).*C1(2,1,:))),1,[]) + C2 - C3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merge the gammas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
w_ =  reshape(w,1, 1,[]);
a_k = sum(w_.*a,3)./wb;

iter = 1;
while iter<20
    iter = iter+1;
    
    h_k = log(a_k)-psi(0,a_k)+c;
    hp_k = 1./a_k - psi(1,a_k);
    hb_k = -1./a_k.^2 - psi(2,a_k);
    
    a_new = a_k-(2.*h_k.*hp_k)./(2.*hp_k.^2-h_k.*hb_k); % Halley's
    
    if abs(a_new-a_k)<1e-2
        a_k = a_new;
        break
    else
        a_k = a_new;
    end
    
    a_k=max(a_k,v_min);
end
a2 = a_k;
b2 = a2./(sum(w_.*a./b,3)./wb);
catch
    a2 = sum(w_.*a,3)./wb;
    b2 = sum(w_.*b,3)./wb;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merge the normals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m2 = sum(m.*repmat(w(:)',n_x,1),2)/wb;

P2 = zeros(n_x,n_x);
for i = 1:N
    P2 = P2 + w(i)*(P(:,:,i)+(m(:,i)-m2)*(m(:,i)-m2)');
end
P2 = P2/wb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merge Inverse Wisharts    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_k = mean(nu,3);

iter = 1;
while iter<100
    iter=iter+1;

    h_k  = d*wb*log(v_k-d-1)...
        -wb*sum(psi(0,(v_k-(d+(1:d)'))./2))...
        +C;
    hp_k = d*wb./(v_k-d-1)...
        -0.5.*wb.*sum(psi(1,(v_k-(d+(1:d)'))./2));
    hb_k = -d*wb./((v_k-d-1).^2)...
        -0.25*wb*sum(psi(2,(v_k-(d+(1:d)'))./2));
   
%     v_new = v_k-h_k/hp_k; % Newtons
    v_new = v_k-(2.*h_k.*hp_k)./(2*hp_k.^2-h_k.*hb_k); % Halley's
    
    if abs(v_new-v_k)<1e-2
        v_k = v_new;
        break
    else
        v_k = v_new;
    end
    v_k=max(v_k,nu_min);
end
nu2=v_k;
V2 = zeros(d,d,length(a2));
for i=1:N
    V2(:,:,:) = V2 + w(i)*(nu_(:,:,:,i)-d-1).*(V(:,:,:,i).\eye(d));
end
% V2 = (V2+V2')/2;
V2 = (V2 + permute(V2, [2 1 3])) / 2;
%%%计算V2的逆
% V2 = (nu2-d-1).*wb.*(V2\eye(d));
% 提取矩阵元素（自动沿第三维度展开）
a = V2(1,1,:);  % 所有切片的(1,1)元素
b = V2(1,2,:);  % 所有切片的(1,2)元素
c = V2(2,1,:);  % 所有切片的(2,1)元素
d_ = V2(2,2,:);  % 所有切片的(2,2)元素

% 计算行列式（逐元素运算）
det_V2 = a.*d_ - b.*c;

% 显式计算逆矩阵（直接应用2×2逆矩阵公式）
inv_V2 = zeros(size(V2));
inv_V2(1,1,:) =  d_ ./ det_V2;  % 新矩阵的(1,1)元素
inv_V2(1,2,:) = -b ./ det_V2;  % 新矩阵的(1,2)元素
inv_V2(2,1,:) = -c ./ det_V2;  % 新矩阵的(2,1)元素
inv_V2(2,2,:) =  a ./ det_V2;  % 新矩阵的(2,2)元素

% V2 = (nu2-d-1)*wb*(V2\eye(d));
% reshape((nu2-d-1)*wb, [1, 1, numel((nu2-d-1)*wb)])
V2 = reshape((nu2-d-1)*wb, [1, 1, numel((nu2-d-1)*wb)]).*inv_V2;

w2 = wb;