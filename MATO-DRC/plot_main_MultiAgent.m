function plot_main_MultiAgent(G_ellip,c_ellip,G_zono,c_zono,Reachset_c,Reachset_G,tau,z0,zN,Solutions,nu,d,W,times,N,Na,points,sol_index,zono_or_ellipsoid)

% zono_or_ellipsoid = 0  for ellipsoid reachable sets
% zono_or_ellipsoid = 1  for zonotope reachable sets

map_fig = plot_map(G_ellip,c_ellip,G_zono,c_zono,z0,zN);

m = 3;
p = (N-1)*m*(d+1);


for a = 1:1:Na
    [tr_figs,Bp_figs] = test_b_to_trajectory_MultiAgent(nu,d,W,times,N,z0{1,a},zN{1,a},Solutions((a-1)*p+1:a*p,sol_index),tau,Reachset_c,Reachset_G,zono_or_ellipsoid);
end

%%% the start and end points of trajectories
% for i = 1:1:Na  
%     p0N = plot3(points(i,1),points(i,2),points(i,3),'square','MarkerSize',7,'MarkerEdgeColor','b','MarkerFaceColor','b'); hold on
% end
% 
% hh=legend(map_fig,'Obstacle');
% set(hh,'Interpreter','latex')
% set(hh,'Location','northwest')
% 
% hh=legend([map_fig, Bp_figs{1}, Bp_figs{2}, Bp_figs{3}, Bp_figs{4}, tr_figs{1}, tr_figs{2}, tr_figs{3}, tr_figs{4}, p0N],'Obstacle', 'BP: 1st', 'BP: 2nd', 'BP: 3rd', 'BP: 4th', 'Tr: 1st', 'Tr: 2nd', 'Tr: 3rd', 'Tr: 4th', 'Start/end');
% set(hh,'Interpreter','latex')
% set(hh,'Location','northwest')

% [c_bar] = color_bar();
% colormap(c_bar);
% cb = colorbar('Ticks',[0, 0.25, 0.5, 0.75, 1],'TickLabels',{'$0s$','$5s$','$10s$','$15s$','$20s$'},'TickLabelInterpreter','latex','Location','east');
% cb.Label.String = 'Reachable sets with time';
% cb.Label.Interpreter = 'latex';


axis([-2 42 -2 42 -2 42]);

% set(gcf,'position',[x0,y0,width,height]);
set(gcf,'position',[50,50,800,600]);
% set(gcf,'position',[50,50,400,300]);

set(gcf,'defaultfigurecolor','w'); 
% set(groot,'defaultLineLineWidth',0.5)

% ti=title('$k=8$');
% set(ti,'Interpreter','latex')
xl=xlabel('$r_1$');
set(xl,'Interpreter','latex')
yl=ylabel('$r_2$');
set(yl,'Interpreter','latex')
zl=zlabel('$r_3$');
set(zl,'Interpreter','latex')
t_size = 12;
set(gca,'FontSize',t_size);

set(gca, 'LooseInset', [0,0,0,0]);  %%%% Remove additional blank borders

% view(az,el); %%% az [0 360] el [-90 90]
view(-32,12);
% view(0,0);  %% r1r3
% view(0,90); %% r1r2

set(gcf, 'PaperPositionMode','auto');

% set(gcf, 'PaperSize', [8.4, 6.3]);
% set(gcf, 'PaperSize', [4.2, 3.15]);

% print(gcf,'-dpdf','-r400','./5Agents_4N_with_BezierPolygon_r1r2');
% print(gcf,'-dpdf','-r400','./5Agents_4N_with_reachableSet');
% print(gcf,'-dpdf','-r400','./5Agents_4N_with_reachableSet_zono_r1r3');
% print(gcf,'-dpdf','-r400','./5Agents_4N_largeGap_2');

end

function map_fig = plot_map(G_ellip,c_ellip,G_zono,c_zono,z0,zN)

figure(1);
color = [143, 209, 225];
color = color /255;
EdgeColor = color;


alf = 0.2;
Edgealf = 1;
LineWidth = 0.3;
LineStyle = '-';

accuracy = 24;

c_cz = c_zono{1};   %%%% zonotope 1
G_cz = G_zono{1};
Ae_cz = [];
be_cz = [];
plot_constrained_zonotope_3D(c_cz,G_cz,Ae_cz,be_cz,color,alf,LineWidth,LineStyle,EdgeColor,Edgealf);

c_cz = c_zono{2};  %%%% zonotope 2
G_cz = G_zono{2};
Ae_cz = [];
be_cz = [];
plot_constrained_zonotope_3D(c_cz,G_cz,Ae_cz,be_cz,color,alf,LineWidth,LineStyle,EdgeColor,Edgealf);

c_ep = c_ellip{1}; %%%% ellipsoid 1
G_ep = G_ellip{1};
plot_ellipsoid_3D(c_ep,G_ep,color,alf,LineWidth,LineStyle,EdgeColor,Edgealf,accuracy);

c_ep = c_ellip{2};  %%%% ellipsoid 2
G_ep = G_ellip{2};
map_fig = plot_ellipsoid_3D(c_ep,G_ep,color,alf,LineWidth,LineStyle,EdgeColor,Edgealf,accuracy);

% h_cy = 36;
% c_cy = [30;30];
% G_cy = [5 0;0 5];
% plot_cylinder_3D(h_cy,c_cy,G_cy,color,alf,LineWidth,LineStyle);

axis([-5 45 -5 45 -5 45]);

% set(fig1,'LineWidth',0.01,'LineStyle',':');

%%%%% start and end point
% plot3(z0(1), z0(4), z0(7), 'square','MarkerFaceColor','b'); hold on
% plot3(zN(1), zN(4), zN(7), 'pentagram','MarkerFaceColor','r');  hold on

end

function plot_trajectory()
end

function plot_reachability_set()
end

function [colort] = color_t(t)
%%%% gradient color
colors = cell(5,1);

% colors{1,1} = [31 120 180] / 255;
% colors{2,1} = [77 174 74] / 255;
% colors{3,1} = [253 141 60] / 255;
% colors{4,1} = [188 128 189] / 255;
% colors{5,1} = [228 26 28] / 255;

% colors{1,1} = [230 240 254] / 255;
% colors{2,1} = [185 210 243] / 255;
% colors{3,1} = [130 170 231] / 255;
% colors{4,1} = [73 108 206] / 255;
% colors{5,1} = [26 49 139] / 255;

% colors{1,1} = [252 187 161] / 255;
% colors{2,1} = [252 146 114] / 255;
% colors{3,1} = [251 106 74] / 255;
% colors{4,1} = [239 59 44] / 255;
% colors{5,1} = [203 24 29] / 255;

% colors{1,1} = [173,191,251] / 255;
% colors{2,1} = [132,161,249] / 255;
% colors{3,1} = [88,130,248] / 255;
% colors{4,1} = [49,99,235] / 255;
% colors{5,1} = [7,72,205] / 255;

% colors{1,1} = [30,144,255] / 255;
% colors{2,1} = [135,206,250] / 255;
% colors{3,1} = [238,130,238] / 255;
% colors{4,1} = [255,105,180] / 255;
% colors{5,1} = [220,20,60] / 255;

colors{1,1} = [30,144,255] / 255;
colors{2,1} = [0,191,255] / 255;
colors{3,1} = [135,206,250] / 255;
colors{4,1} = [255,105,180] / 255;
colors{5,1} = [220,20,60] / 255;



% times = [0;5;10;15;20];
i = fix(t / 5);
b = (t - i * 5) / 5;
if i == 4
    colort = colors{5,1};
else
    colort = colors{i+1,1} + (colors{i+2} - colors{i+1}) * b;
end
end

function [c_bar] = color_bar()
%%%% gradient color
dt = 0.5;
t0 = 0;
tN = 20;
ni = tN / dt;
c_bar = zeros(ni+1,3); % 
for i = 0:1:ni
    t = t0 + i * dt;
    c_bar(i+1,:) = color_t(t);
end
end