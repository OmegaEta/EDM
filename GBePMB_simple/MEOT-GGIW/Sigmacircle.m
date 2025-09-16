function [x,y,H] = Sigmacircle(cx,cy,P,N,color,style,isfill)
%cx cy��Բ�������꣬P��Բ��״��N����׼�color style ��ͼ����
if nargin < 4 %��������ĸ���
    N=1;
    color = 2;
    style ='-'; 
elseif nargin < 5
    color = 2;
    style ='-';
elseif nargin < 6
    style ='-';
end
N2 = 20;
H = [];
num = size(cx,1);
for i = 1 : num
    if size(P(:,:,i),1) == 2
        sqrtP = N*sqrtm_2by2(P(:,:,i));
    else
        sqrtP = N*sqrtm(P(:,:,i));
    end

    phi = 0:pi/N2:2*pi;
    xy = sqrtP*[cos(phi); sin(phi)];
    xy = xy + [cx(i)*ones(1,length(phi)) ; cy(i)*ones(1,length(phi))];

    x = xy(1,:)';
    y = xy(2,:)';
    
    c = ['r' 'g' 'b' 'c' 'm' 'y' 'k' 'w'];
    %    ��   ��  �� ���� ��� ��  ��  ��
    hp = plot(x,y,strcat(style,c(mod(color-1,7)+1)),'linewidth',1);
    H = [H;hp];

    if isfill ==1 
        h = fill(x,y,c(mod(color-1,7)+1));
        set(h,'facealpha',0.1);
        H = [H;h];
    end
    hold on;
%     pbaspect([1 1 1]);
%     axis([-200 200 -200 200])
end
end