function [Z,groundTruth] = CreateMeasFunction(i,i0,D,model)
% 
if i0==0
    if i==1 %ƽ��
        Scenario_2;
    elseif i==2 %����
        Scenario_3;
    else  %28��Ŀ��
        Scenario_5;
    end
else
    if i==1
        
    elseif i==2
        Z = load('./measures/Z2.mat').Z;
        groundTruth = load('./measures/groundTruth2.mat').groundTruth;
    else
        Scenario_4;
    end
end
end

