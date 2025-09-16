
m = 100;
GOSPA_imp = zeros([100 4 m]);
GOSPA_org = zeros([100 4 m]);
n = 0;
for i = 1:6
    filename1 = strcat('./measures/新气象2/GOSPA_impr3_',string(i),'.mat');
    filename2 = strcat('./measures/新气象2/GOSPA_orig3_',string(i),'.mat');
    G1 = load(filename1);
    G2 = load(filename2);
    n=n+1;
    GOSPA_imp(:,:,n) = G1.GOSPA_impr;
    GOSPA_org(:,:,n) = G2.GOSPA_orig;
end

for i = 1:5
    filename1 = strcat('./measures/新气象/GOSPA_impr3_',string(i),'.mat');
    filename2 = strcat('./measures/新气象/GOSPA_orig3_',string(i),'.mat');
    G1 = load(filename1);
    G2 = load(filename2);
    n=n+1;
    GOSPA_imp(:,:,n) = G1.GOSPA_impr;
    GOSPA_org(:,:,n) = G2.GOSPA_orig;
end

for i = 1:8
    filename1 = strcat('./measures/新气象3/GOSPA_impr3_',string(i),'.mat');
    filename2 = strcat('./measures/新气象3/GOSPA_orig3_',string(i),'.mat');
    G1 = load(filename1);
    G2 = load(filename2);
    n=n+1;
    GOSPA_imp(:,:,n) = G1.GOSPA_impr;
    GOSPA_org(:,:,n) = G2.GOSPA_orig;
end

for i = 1:12
    filename1 = strcat('./measures/新气象4/GOSPA_impr3_',string(i),'.mat');
    filename2 = strcat('./measures/新气象4/GOSPA_orig3_',string(i),'.mat');
    G1 = load(filename1);
    G2 = load(filename2);
    n=n+1;
    GOSPA_imp(:,:,n) = G1.GOSPA_impr;
    GOSPA_org(:,:,n) = G2.GOSPA_orig;
end
for i = 1:25
    filename1 = strcat('./measures/新气象5/GOSPA_impr3_',string(i),'.mat');
    filename2 = strcat('./measures/新气象5/GOSPA_orig3_',string(i),'.mat');
    G1 = load(filename1);
    G2 = load(filename2);
    n=n+1;
    GOSPA_imp(:,:,n) = G1.GOSPA_impr;
    GOSPA_org(:,:,n) = G2.GOSPA_orig;
end


for i = 1:8
    filename1 = strcat('./measures/新气象3/GOSPA_impr3_',string(i),'.mat');
    filename2 = strcat('./measures/新气象3/GOSPA_orig3_',string(i),'.mat');
    G1 = load(filename1);
    G2 = load(filename2);
    n=n+1;
    GOSPA_imp(:,:,n) = G1.GOSPA_impr;
    GOSPA_org(:,:,n) = G2.GOSPA_orig;
end

GOSPA02_improved = sum(GOSPA_imp,3) / n;
GOSPA02_original = sum(GOSPA_org,3) / n;
figure(6);
clf(6);
c = 1;
plot(1:100,GOSPA02_improved(:,c),'r','LineWidth',1.5);
hold on;
plot(1:100,GOSPA02_original(:,c),'b','LineWidth',1.5);
legend('改进PMBM','原PMBM');
xlabel('time');
ylabel('GOSPA');

figure(7);
clf(7);
c = 2;
plot(1:100,GOSPA02_improved(:,c),'r','LineWidth',1.5);
hold on;
plot(1:100,GOSPA02_original(:,c),'b','LineWidth',1.5);
legend('改进PMBM','原PMBM');
xlabel('time');
ylabel('location error');

Scenario_5;
N_t = zeros(100,1);
for t=1:100
    N_t(t) = size(groundTruth{t}.x,2);
end

N_improved= N_t-GOSPA02_improved(:,3)+GOSPA02_improved(:,4);
N_original = N_t-GOSPA02_original(:,3)+GOSPA02_original(:,4);
figure(8);
clf(8);
hold on;
box on;
plot(1:100,N_improved,'r','LineWidth',1.5);
plot(1:100,N_original,'b','LineWidth',1.5);
plot(1:100,N_t,'k:','LineWidth',1.5);
xlabel('time(s)');
ylabel('已知目标集的基数');
legend('改进PMBM','原始PMBM','真实值');