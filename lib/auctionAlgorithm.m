function [Customer2Item,Item2Customer]=auctionAlgorithm(RewardMatrix)
%拍卖算法实现分配最大收益
epsilon=0.1; %amount of deviation from the optimal reward

[NofCustomers,NofItems]=size(RewardMatrix);
if (NofCustomers>NofItems)
    error('Number of columns must be greater than or equal to the number of rows');
end

Item2Customer=zeros(1,NofItems);
Customer2Item=zeros(1,NofCustomers);

% 初始化计数器，防止无限循环
max_iterations = 1000;
iteration_count = 0;

while ~isempty(find(Customer2Item==0, 1)) && iteration_count < max_iterations
    iteration_count = iteration_count + 1;
    
    assigned = false; % 标记本次循环是否有新分配
    
    if (NofItems==1) % 如果只有一个物品
        [~, Item2Customer]=max(RewardMatrix); % 分配物品给最佳客户
        Customer2Item(Item2Customer)=1; % 标记客户已分配
        assigned = true;
    else
        for i=1:NofCustomers
            if ~Customer2Item(i) % 如果客户未分配
                [maxval,maxind]=max(RewardMatrix(i,:)); % 找到最大收益物品
                RewardMatrix(i,maxind)=min(RewardMatrix(i,:))-1; % 临时降低最大收益，以便找到次大值
                [secondmaxval,secondmaxind]=max(RewardMatrix(i,:)); % 找到次大收益物品
                RewardMatrix(i,maxind)=maxval; % 恢复最大收益
                
                if maxval > -inf % 确保存在可分配的物品
                    Customer2Item(i)=maxind; % 分配物品给客户
                    if Item2Customer(maxind) % 如果物品已被分配
                        Customer2Item(Item2Customer(maxind))=0; % 取消原客户的分配
                    end
                    Item2Customer(maxind)=i; % 更新物品分配
                    
                    % 降低物品对其他客户的吸引力
                    if maxval > secondmaxval % 防止负值
                        RewardMatrix(:,maxind)=RewardMatrix(:,maxind)-(maxval-secondmaxval+epsilon);
                    end
                    
                    assigned = true;
                end
            end
        end
    end
    
    % 如果本次循环没有任何分配，说明无法继续分配，退出循环
    if ~assigned
        break;
    end
end

% 检查是否达到最大迭代次数
if iteration_count >= max_iterations
    warning('达到最大迭代次数，可能存在分配问题');
    keyboard
end