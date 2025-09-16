clear;clc
dbstop if error
addpath(genpath('..\Scenario'));addpath(genpath('..\lib'));addpath(genpath('..\lib_EDM'))

scenario = 1;
ii=0;
%特点场景中特点目标在特点时刻的间距，即一号场景中，两个平行目标间距，预设5。
%二号场景中，两个时刻50时，初始两目标的间距，预设12
if scenario == 1
    modelparas1;
else
    D = 12;
end
%Parameter setting
modelparas1;
%Number of Monte Carlo Simulations
numMC = 100;
%Parameters used in GOSPA metric
c = 20;
p = 1;
%Number of time steps
K = model.K;

GOSPA_improved = zeros(K,4,numMC);
GOSPA_original = zeros(K,4,numMC);
trajectoryEstimates_improved = cell(numMC,1);
trajectoryEstimates_original = cell(numMC,1);
simulation_time = zeros(numMC,1);

for t = 1:numMC
    % Initialisation
    PPP_improved.w = log(model.birth.w);
    PPP_improved.GGIW = model.birth.GGIW;
    PPP_original.w = log(model.birth.w);
    PPP_original.GGIW = model.birth.GGIW;
    
    MBM_improved.w = [];     % Global hypotheses weights
    MBM_improved.track = {}; % Locl hypotheses trees
    MBM_improved.table = []; % Global hypotheses look-up table
    MBM_original.w = [];     % Global hypotheses weights
    MBM_original.track = {}; % Locl hypotheses trees
    MBM_original.table = [];
    
    estimates_improved = cell(K,1);
    trajectoryEstimates_improved{t} = cell(K,1);
    estimates_original = cell(K,1);
    trajectoryEstimates_original{t} = cell(K,1);
    
    %即时产生量测
%     [Z,groundTruth] = CreateMeasFunction(scenario,ii,D,model);
    
    cost_time = zeros(K,1);
    tic
    for k = 1:K
        [t,k]
        tic;
        %Update step
        [PPP_improved,MBM_improved] = updatePMBM_improved(PPP_improved,MBM_improved,Z{k},k,model);%改进算法
%         [PPP_original,MBM_original] = updatePMBM_original(PPP_original,MBM_original,Z{k},k,model);
        %Extract estimates (both estimate of the current time and the
        %estimate of the full trajectory) 
        [estimates_improved{k},trajectoryEstimates_improved{t}{k}] = estimatorks(MBM_improved,model);
        [estimates_original{k},trajectoryEstimates_original{t}{k}] = estimatorks(MBM_original,model);
        %Evaluate filtering performance using GOSPA
        GOSPA_improved(k,:,t) = GOSPAmetric_ks(estimates_improved{k},groundTruth{k},c,p);
        GOSPA_original(k,:,t) = GOSPAmetric_ks(estimates_original{k},groundTruth{k},c,p);
        
        % 画出跟踪结果
       if 1 %%t==0 && k>=40 && k<=55&& mod(k,2)==0
             figure(2);
             box on;
%               clf(2);
             axis([-200,200,-200,200]);
             hold on;
             plot(Z{k}(1,:),Z{k}(2,:),"k.");
             i0=size(estimates_improved{k}.g,2);
             for j = 1:i0
                 figure(2);
                
                 [x, y] = Sigmacircle(estimates_improved{k}.x(1,j),estimates_improved{k}.x(2,j),estimates_improved{k}.X(:,:,j),2,1,'-',1);
             end
             drawnow;
             i0=size(estimates_original{k}.g,2);
             for j = 1:i0
                 figure(2);
%                  [x, y] = Sigmacircle(estimates_original{k}.x(1,j),estimates_original{k}.x(2,j),estimates_original{k}.X(:,:,j),2,3,'--',1);
             end
             drawnow;
        end

         %Prediction Step
        if k < K
            [PPP_improved,MBM_improved] = predictPMBM(PPP_improved,MBM_improved,model);
        end
        if k < K
%             [PPP_original,MBM_original] = predictPMBM(PPP_original,MBM_original,model);
        end
        cost_time(k) = toc;
    end
    GOSPA_impr = GOSPA_improved(:,:,t);
    GOSPA_orig = GOSPA_original(:,:,t);
    
    filename1 = strcat("./measures/GOSPA_impr",string(scenario),"_",string(t+25),".mat");
    filename2 = strcat("./measures/GOSPA_orig",string(scenario),"_",string(t+25),".mat");
    save(filename1,'GOSPA_impr');
    save(filename2,'GOSPA_orig');
    
    filename3 = strcat("./MBM/trajectoryEstimates_improved",string(scenario),"_",string(t+25),".mat");
    filename4 = strcat("./MBM/trajectoryEstimates_original",string(scenario),"_",string(t+25),".mat");
    trajectory_improved = trajectoryEstimates_improved{t};
    trajectory_original = trajectoryEstimates_original{t};
    save(filename3,'trajectory_improved');
    save(filename4,'trajectory_original');
    
     simulation_time(t) = toc;
    
      %每十个时刻画出GOSPA的平均
    if 1 %mod(t,5)==1
        GOSPA02_improved = sum(GOSPA_improved,3) / t;
        GOSPA02_original = sum(GOSPA_original,3) / t;
        figure(4);
        clf(4);
        plot(1:100,GOSPA02_improved(:,1),'r');
        hold on;
        plot(1:100,GOSPA02_original(:,1),'b');
        legend('改进PMBM','原PMBM');
        xlabel('time');
        ylabel('GOSPA');
    end
end
 %100个时刻的GOSPA的平均
GOSPA02_improved = sum(GOSPA_improved,3) / numMC;
GOSPA02_original = sum(GOSPA_original,3) / numMC;