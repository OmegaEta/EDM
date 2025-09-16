function drawMeas(XY,color,style,fig,clf_)
%clarify 
%clf_ 有值为清空图；无值为保留原图数据
if nargin < 5 %输入参数的个数
    hold on;
elseif nargin < 6
    figure(fig);
    clf(fig);
end

c = ['r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k'];
    %红   绿  蓝 青绿 洋红 黄  黑 红   绿  蓝 青绿 洋红 黄  黑 
N = size(XY,1);
for i = 1:N
    figure(fig);
%     clf(fig);
    pause(0.5);
    if color == 0
        cc = i; 
    else
        cc = color; 
    end
    if iscell(XY)
        if isempty(XY{i})
            continue;
        end
        plot(XY{i}(:,1),XY{i}(:,2),strcat(style,c(mod(cc-1,7)+1)));
    else
%         plot(XY(i,1),XY(i,2),strcat(style,c(mod(cc-1,7)+1)));
        plot(XY(:,1),XY(:,2),strcat(style,c(mod(cc-1,7)+1)));
        break;
    end
%     hold on;
end
    
end

