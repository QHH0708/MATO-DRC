function plot_main(G_ellip,c_ellip,G_zono,c_zono,z0,zN,Solutions,nu,d,W,times,N,sol_index)

maplegend = plot_map(G_ellip,c_ellip,G_zono,c_zono);

tra_color = [213, 33, 32] / 255;  %%%% red
BP_color = [252, 158, 127] / 255;
BP_Edgecolor = [0.4940 0.1840 0.5560];

num_s = size(sol_index,2);
for i = 1:1:num_s
    [trajectoryfig,Bpolygonfig] = plot_traj_Bezier_polygon(nu,d,W,times,N,z0,zN,Solutions(:,sol_index(i)),tra_color,BP_color,BP_Edgecolor);
end

%%%% start and end point
% plot3(z0(1), z0(4), z0(7), 'square','MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor','b'); hold on
% plot3(zN(1), zN(4), zN(7), 'pentagram','MarkerSize',8,'MarkerEdgeColor','r','MarkerFaceColor','r');  hold on

p0 = plot3(z0(1), z0(4), z0(7), 'square','MarkerSize',7,'MarkerEdgeColor','b','MarkerFaceColor','b'); hold on
pN = plot3(zN(1), zN(4), zN(7), 'pentagram','MarkerSize',7,'MarkerEdgeColor','r','MarkerFaceColor','r');  hold on

axis([-2 42 -2 42 -2 42]);

% hh=legend([maplegend trajectoryfig p0 pN],'Obstacle','Trajectory','Start','End');
hh=legend([maplegend Bpolygonfig trajectoryfig p0 pN],'Obstacle','Bezier polygon','Trajectory','Start','End');
set(hh,'Interpreter','latex')
set(hh,'Location','northwest')

% set(gcf,'position',[x0,y0,width,height]);
set(gcf,'position',[50,50,800,600]);

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

% view(az,el); %%% az [0 360] el [-90 90]
view(-39,11);

% print(gcf,'-dpdf','-r600','./1UAV_with_Bezier');

end

function [maplegend] = plot_map(G_ellip,c_ellip,G_zono,c_zono)
figure(1);

color = [143, 209, 225];
color = color /255;
EdgeColor = color;


alf = 0.2;
Edgealf = 0.8;
LineWidth = 0.2;
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
maplegend = plot_ellipsoid_3D(c_ep,G_ep,color,alf,LineWidth,LineStyle,EdgeColor,Edgealf,accuracy);

end
