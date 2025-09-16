function polt_PDstd(PD,c,if_reverse)

std1 = std(PD);

figure;clf;hold on;
yyaxis left;
set(gca,'ytick',[],'yticklabel',[]);

yyaxis right;
stem(std1);
colororder(c);
set(gca,'xtick',[],'xticklabel',[])

if if_reverse
    ax = gca;
    ax.YAxis(2).Direction = 'reverse';
end

end