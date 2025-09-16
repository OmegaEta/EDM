clear;clc;


s = 1;
T = 100;
N = 90;

NUM_improved = zeros(T,N);
NUM_original = zeros(T,N);
NUM_A=zeros(T,N);
for n = 1:N
    filename1 = strcat('./MBM/ƽ�г����켣/trajectoryEstimates_improved',string(s),'_',string(n),'.mat');
    filename2 = strcat('./MBM/ƽ�г����켣/trajectoryEstimates_original',string(s),'_',string(n),'.mat');
    load(filename1)
    load(filename2);
    for i=1:100
        NUM_A(i,n)=2;
        if(isempty(trajectory_improved{i}))
            NUM_improved(i,n) = 0;
        else
            idx_improved = [[trajectory_improved{i}.t_death]>=i & [trajectory_improved{i}.t_birth]<=i];
            NUM_improved(i,n) = sum(idx_improved);
        end
        if(isempty(trajectory_original{i}))
            NUM_original(i,n) = 0;
        else
            idx_original = [[trajectory_original{i}.t_death]>=i & [trajectory_original{i}.t_birth]<=i];
            NUM_original(i,n) = sum(idx_original);
        end
    end
end
ave_A = sum(NUM_A,2)/N;
ave_improved = sum(NUM_improved,2)/N;
ave_original = sum(NUM_original,2)/N;
figure(2);
hold on;
box on;
plot([1:100],ave_improved,'r','LineWidth',1.5);
plot([1:100],ave_original,'b','LineWidth',1.5);
plot([1:100],ave_A,'black','LineWidth',1.5);
xlabel('time(s)');
ylabel(['Object Cardinality']);
legend('KS-GGIW-PMBM','GGIW-PMBM','Real Value');