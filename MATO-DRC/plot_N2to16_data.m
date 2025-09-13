

N = 2:2:16;

%%%% computational_time / solution_number
CT = [ct_n2/sn_n2, ct_n4/sn_n4, ct_n6/sn_n6, ct_n8/sn_n8, ct_n10/sn_n10, ct_n12/sn_n12, ct_n14/sn_n14, ct_n16/sn_n16];
%%%% objective value
OV = [obj_val_n2, obj_val_n4, obj_val_n6, obj_val_n8, obj_val_n10, obj_val_n12, obj_val_n14, obj_val_n16];


figure;
% set(gcf,'defaultfigurecolor','w');

yyaxis left; 
plot(N, CT, 'o--','color', 'b', 'MarkerSize', 8, 'LineWidth', 2);
yl=ylabel('Computational time/s');
set(yl,'Interpreter','latex')
ylim([0, 30]);

yyaxis right;     
plot(N, OV, 's-','color', 'r', 'MarkerSize', 8,  'LineWidth', 2);
yl=ylabel('Objective value');
set(yl,'Interpreter','latex')
ylim([1.5, 3]);    


xl=xlabel('$N$');
set(xl,'Interpreter','latex')
xlabel('N', 'FontSize', 12);
% xlim([-0.01, 16.01]);
grid on;
hh = legend('Computational Time', 'Objective Value');
set(hh,'Interpreter','latex')
set(hh,'Location','northwest')

set(gca, 'Layer', 'bottom');

t_size = 12;
set(gca,'FontSize',t_size);
set(gcf,'position',[50,50,800,400]);

set(gcf, 'PaperPositionMode','auto');
set(gcf, 'PaperSize', [8.4, 6.3]);
% print(gcf,'-dpdf','-r400','./comptime_objval');



