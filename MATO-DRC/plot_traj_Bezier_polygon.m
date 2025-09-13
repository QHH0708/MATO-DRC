function [trajectoryfig,Bpolygonfig] = plot_traj_Bezier_polygon(nu,d,W,times,N,z0,zN,zd,tra_color,BP_color,BP_Edgecolor)
%%% zd is a vector with (N-1)*m*(d+1) elements
%%% zd denotes the parameters of interior-point constraints

sum_nu = sum(nu);
m = size(nu,1);
bar_nu_z = 2*max(nu);
[bold_M_new,zd_matrix_bold_b,vector_z0N] = bold_M_matrix(nu,d,W,times,N,z0,zN);
zd_new = zd;
bold_b_new = zd_matrix_bold_b * zd_new + vector_z0N;  %%% only applicable to d = 0 case
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
% color = [213, 33, 32] / 255;  %%%% red
% color = [155 34 39]/255;
% color = [0 159 105]/ 255; %%%% green
% color = [77 174 74] / 255; %%% green
% color = 'g';
color = tra_color;
trajectoryfig = plot3(Points(1,:),Points(2,:),Points(3,:),'LineWidth',1,'Color',color); hold on
end


function [ret] = polynomial_vector(order,t)
%%%%%% order = the degree of the trajectory + 1
ret = zeros(1, order);
for i = 1:1:order
    ret(1,i) = t^(i - 1) / factorial(i - 1);
end
end