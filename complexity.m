%complexity analysis
clear;
close all;
M_n = 5;
comp = zeros(M_n, 3);
Ms = [9;  49; 100;  225; 400];
for m = 1 : M_n
    comp(m, 1) = Ms(m, 1) ^ 3;
    comp(m, 2) = comp(m, 1) * 1.21 * 1e4;
    comp(m, 3) = comp(m, 1) * 1.4641 * 1e8;
end
comp = comp + 5e6;
comp = log10(comp);
h = figure;
set(h,'PaperType','A4');
axes('FontSize',16);
bar(1:M_n, comp)
le = legend('Proposed', 'DISPARE [37]','Subspace [40]', 'Location','Northeast');
set(le,'Fontsize',14,'Fontname','Times')
set(gca,'XTickLabel',{Ms})
ylim([0, max(max(comp))+5])
grid on
xlabel('Number of the BS antennas \it M','Fontsize',16,'Fontname','Times')
ylabel('Base 10 logarithm of the complexity','Fontsize',16,'Fontname','Times')
% print(h,'-dpdf','complexity_antenna')