
clear;clc;

drawN;
drawLocationerror;

ave_i = ave_improved1./ave_improved;
ave_o = ave_original1./ave_original;

figure(3);
hold on;
box on;
plot([1:100],ave_i,'r','LineWidth',1.5);
plot([1:100],ave_o,'b','LineWidth',1.5);
xlabel('time(s)');
ylabel('average location error');
legend('改进PMBM','原始PMBM');