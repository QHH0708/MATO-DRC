function plot_main_dif_N(G_ellip,c_ellip,G_zono,c_zono,z0,zN,nu,d,W,sN2,iN2,sN4,iN4,sN8,iN8,sN16,iN16)
%%%%% s: Solutios,  i：sol_index


maplegend = plot_map(G_ellip,c_ellip,G_zono,c_zono);


t0 = 0;
tN = 20;


%%%% tra_color: trajectory color,  BP_color: Bezier polygon color
N=2;
times = zeros(N+1,1);
for n = 1:1:N
    times(n+1) = (tN/N) * n;
end
tra_color = [49, 124, 183] / 255;  
BP_color = [182, 215, 232] / 255;
BP_Edgecolor = [109, 173, 209] / 255;
for i = 1:1:iN2
    [tfigN2,BfigN2] = test_b_to_trajectory_color(nu,d,W,times,N,z0,zN,sN2(:,i),tra_color,BP_color,BP_Edgecolor);
end

N=4;
times = zeros(N+1,1);
for n = 1:1:N
    times(n+1) = (tN/N) * n;
end
tra_color = [213, 33, 32] / 255;  %%%% red
BP_color = [252, 158, 127] / 255;
BP_Edgecolor = [0.4940 0.1840 0.5560];
for i = 1:1:iN4
    [tfigN4,BfigN4] = test_b_to_trajectory_color(nu,d,W,times,N,z0,zN,sN4(:,i),tra_color,BP_color,BP_Edgecolor);
end

N=8;
times = zeros(N+1,1);
for n = 1:1:N
    times(n+1) = (tN/N) * n;
end
tra_color = [145, 60, 137] / 255;  
% BP_color = [245, 232, 238] / 255;
BP_color = [246, 103, 160] / 255;
BP_Edgecolor = [194, 144, 197] / 255;
for i = 1:1:iN8
    [tfigN8,BfigN8] = test_b_to_trajectory_color(nu,d,W,times,N,z0,zN,sN8(:,i),tra_color,BP_color,BP_Edgecolor);
end

N=16;
times = zeros(N+1,1);
for n = 1:1:N
    times(n+1) = (tN/N) * n;
end
tra_color = [27, 124, 61] / 255;  
BP_color = [178, 225, 185] / 255;
BP_Edgecolor = [119, 202, 197] / 255;
for i = 1:1:iN16
    [tfigN16,BfigN16] = test_b_to_trajectory_color(nu,d,W,times,N,z0,zN,sN16(:,i),tra_color,BP_color,BP_Edgecolor);
end

%%%% 起点和终点
p0 = plot3(z0(1), z0(4), z0(7), 'square','MarkerSize',7,'MarkerEdgeColor','b','MarkerFaceColor','b'); hold on
pN = plot3(zN(1), zN(4), zN(7), 'pentagram','MarkerSize',7,'MarkerEdgeColor','r','MarkerFaceColor','r');  hold on

axis([-3 42 -3 42 -3 42]);

% hh=legend([maplegend trajectoryfig p0 pN],'Obstacle','Trajectory','Start','End');
hh=legend([maplegend BfigN2 BfigN4 BfigN8 BfigN16 tfigN2 tfigN4 tfigN8 tfigN16 p0 pN],'Obstacle','Bp: N=2', 'Bp: N=4', 'Bp: N=8', 'Bp: N=16', 'Tr: N=2','Tr: N=4','Tr: N=8','Tr: N=16','Start','End');
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

set(gca, 'LooseInset', [0,0,0,0]);  %%%% 去除额外的空白边框

set(gcf, 'PaperPositionMode','auto');

set(gcf, 'PaperSize', [8.4, 6.3]);

% view(az,el); %%% az [0 360] el [-90 90]
view(-32,12);

print(gcf,'-dpdf','-r400','./1Agent_N24816_Bezier');

end

function [maplegend] = plot_map(G_ellip,c_ellip,G_zono,c_zono)
figure(1);
% color = [154, 197, 244]; %% b
% color = [27, 120, 178]; %% b
% color = [34, 189, 210];
% color = [134, 148 173];
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
maplegend = plot_ellipsoid_3D(c_ep,G_ep,color,alf,LineWidth,LineStyle,EdgeColor,Edgealf,accuracy);

end


function [trajectoryfig,Bpolygonfig] = test_b_to_trajectory_color(nu,d,W,times,N,z0,zN,zd,tra_color,BP_color,BP_Edgecolor)
%%% zd (N-1) 行  m*(d+1) 列 矩阵 估计中间点约束
sum_nu = sum(nu);
m = size(nu,1);
bar_nu_z = 2*max(nu);
[bold_M_new,zd_matrix_bold_b,vector_z0N] = bold_M_matrix(nu,d,W,times,N,z0,zN);
zd_new = zd;
bold_b_new = zd_matrix_bold_b * zd_new + vector_z0N;  %%% 只适用于 d = 0 的情况
bold_c_all = bold_M_new \ bold_b_new;

select_cn = cell(N,1);
for i = 1:1:N
    select_cn_temp = zeros(2*sum_nu,2*N*sum_nu);
    select_cn_temp(:,(i-1)*2*sum_nu + 1:i*2*sum_nu) = eye(2*sum_nu);
    select_cn{i,1} = select_cn_temp;
end

%%%%% value of objective function
Wint = Wint_matrix_ofobjectivefunction(nu,W,times,N);
% disp(bold_c_all);
val_objfunction = bold_c_all' * Wint * bold_c_all;
disp(val_objfunction);
[reorganize_matrices] = reorganize_each_flatoutput(nu,W);
QB = QB_expanded_matrices_Bernstein(bar_nu_z);
TsQB = cell(N,1); 
for n = 1:1:N
    Ts = Ts_time_scaling_matrix_Bernstein(bar_nu_z,times(n,1),times(n+1,1));
    TsQB{n,1} = Ts' / QB;
end
%%% control_P = zd' * A_control_P + B_control_P
A_control_P = cell(N,m);
c_control_P = cell(N,m); 
for n = 1:1:N
    for i = 1:1:m
        A_control_P{n,i} = zd_matrix_bold_b' / bold_M_new' * select_cn{n}' * reorganize_matrices{i}' * TsQB{n};
        c_control_P{n,i} = vector_z0N' / bold_M_new' * select_cn{n}' * reorganize_matrices{i}' * TsQB{n};
    end
end
for n = 1:1:N   %%% trajectory piece
    bold_c_n = select_cn{n,1} * bold_c_all;
    new_parameters = zeros(m,2*max(nu));
    for i = 1:1:m  %%% for each output / each dimension
        new_parameters(i,:) = bold_c_n' * reorganize_matrices{i,1}';
    end
    trajectoryfig = plot_piece_trajectory(new_parameters,nu,times(n),times(n+1),tra_color);
    Bpolygonfig = 0;
    Bpolygonfig = plot_piece_Bezier_polygon(new_parameters,bar_nu_z,times(n),times(n+1),BP_color,BP_Edgecolor);
    
    A_CP = cell(1,m);
    c_CP = cell(1,m);
    for i = 1:1:m
        A_CP{i} = A_control_P{n,i};
        c_CP{i} = c_control_P{n,i};
    end
end

end

function [B_p] = plot_piece_Bezier_polygon(new_parameters,bar_nu_z,t_start,t_end,BP_color,BP_Edgecolor)
[Ts] = Ts_time_scaling_matrix_Bernstein(bar_nu_z,t_start,t_end);
[QB] = QB_expanded_matrices_Bernstein(bar_nu_z);
control_P = new_parameters * Ts' / QB;

color = BP_color;

alf = 0.4;
EdgeColor = BP_Edgecolor;
Edgealf = 0.8;
LineWidth = 0.2;
LineStyle = '-';

Bezier_polygon = Polyhedron(control_P');
B_p = plot(Bezier_polygon,'color',color); hold on

alpha(B_p,alf);
set(B_p,'LineWidth',LineWidth,'LineStyle',LineStyle,'EdgeColor',EdgeColor,'EdgeAlpha',Edgealf);
end


function [trajectoryfig] = plot_piece_trajectory(new_parameters,nu,t_start,t_end,tra_color)
%%% new_parameters is a matrix, each raw is a parameter for a dimension
%%% nus is a vector, maximal order + 1 of each dimension
%%% t_start and t_end should be the integral multiples of 0.1
dim = size(new_parameters,1);
bar_nu_z = 2 * max(nu);
% num_p = (t_end - t_start) / 0.05 + 1;
num_p = 50;
dt = (t_end - t_start) / num_p;
Points = zeros(dim,num_p);
point = zeros(dim,1);
for i = 1:1:num_p+1
    t = t_start + (i-1)*dt;  
    for j = 1:1:dim
        point(j,1) = new_parameters(j,:) * polynomial_vector(bar_nu_z,t)';
    end
    Points(:,i) = point;
end
%%% dim = 3
% color = [0.6350 0.0780 0.1840];
% color = [213, 33, 32] / 255;  %%%% 红
% color = [155 34 39]/255;
% color = [0 159 105]/ 255; %%%% 绿
% color = [77 174 74] / 255; %%% 绿
% color = 'g';
color = tra_color;
trajectoryfig = plot3(Points(1,:),Points(2,:),Points(3,:),'LineWidth',1,'Color',color); hold on
end


function [ret] = polynomial_vector(order,t)
%%%%%% order = 多项式的最高次 + 1，0 就是常数
ret = zeros(1, order);
for i = 1:1:order
    ret(1,i) = t^(i - 1) / factorial(i - 1);
end
end
