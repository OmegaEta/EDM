% clear;clc;

s = 1;
T = 100;
N = 90;

NUM_improved = zeros(T,N);
NUM_original = zeros(T,N);
for n = 1:N
    filename1 = strcat('./measures/平行场景GOSPA/GOSPA_impr',string(s),'_',string(n),'.mat');
    filename2 = strcat('./measures/平行场景GOSPA/GOSPA_orig',string(s),'_',string(n),'.mat');
    load(filename1);
    load(filename2);
    NUM_improved(:,n) = GOSPA_impr(:,2);
    NUM_original(:,n) = GOSPA_orig(:,2);
    
end

ave_improved1 = sum(NUM_improved,2)/N;
ave_original1 = sum(NUM_original,2)/N;
figure(2);
clf(2);
hold on;
box on;
plot([1:100],ave_improved1,'r','LineWidth',1.5);
plot([1:100],ave_original1,'b','LineWidth',1.5);
xlabel('time(s)');
ylabel('Location Error');
legend('KS-GGIW-PMBM','GGIW-PMBM');