function [trajectoryfig,Bpolygonfig] = test_b_to_trajectory(nu,d,W,times,N,z0,zN,zd)
%%% zd is a vector with (N-1)*m*(d+1) elements
%%% zd denotes the parameters of interior-point constraints

sum_nu = sum(nu);
m = size(nu,1);
bar_nu_z = 2*max(nu);
%%%% compute the matirx bold_M and bold_b
% [bold_M,bold_b] = bold_M_b_matrix(nu,d,W,times,N,z0,zN,zd);
% bold_c_all = bold_M \ bold_b;


[bold_M_new,zd_matrix_bold_b,vector_z0N] = bold_M_matrix(nu,d,W,times,N,z0,zN);
% zd_new = zeros((N-1)*m*(d+1),1);
% for i = 1:1:N-1
%     zd_new((i-1)*m*(d+1)+1:i*m*(d+1),1) = zd(i,:)';
% end
zd_new = zd;
bold_b_new = zd_matrix_bold_b * zd_new + vector_z0N;  %%% only applicable to d = 0 case
bold_c_all = bold_M_new \ bold_b_new;

select_cn = cell(N,1);
for i = 1:1:N
    select_cn_temp = zeros(2*sum_nu,2*N*sum_nu);
    select_cn_temp(:,(i-1)*2*sum_nu + 1:i*2*sum_nu) = eye(2*sum_nu);
    select_cn{i,1} = select_cn_temp;
    % disp(select_cn{i,1});
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
    % bold_c_n = bold_c_all((n-1)*2*sum_nu + 1:n*2*sum_nu,1);
    bold_c_n = select_cn{n,1} * bold_c_all;
    new_parameters = zeros(m,2*max(nu));
    for i = 1:1:m  %%% for each output / each dimension
        new_parameters(i,:) = bold_c_n' * reorganize_matrices{i,1}';
    end
    % disp(new_parameters);
    trajectoryfig = plot_piece_trajectory(new_parameters,nu,times(n),times(n+1));
    Bpolygonfig = 0;
    Bpolygonfig = plot_piece_Bezier_polygon(new_parameters,bar_nu_z,times(n),times(n+1));
    
    A_CP = cell(1,m);
    c_CP = cell(1,m);
    for i = 1:1:m
        A_CP{i} = A_control_P{n,i};
        c_CP{i} = c_control_P{n,i};
    end
end

end

function [B_p] = plot_piece_Bezier_polygon(new_parameters,bar_nu_z,t_start,t_end)
[Ts] = Ts_time_scaling_matrix_Bernstein(bar_nu_z,t_start,t_end);
[QB] = QB_expanded_matrices_Bernstein(bar_nu_z);
control_P = new_parameters * Ts' / QB;
% disp(control_P)
% control_P = [control_P control_P(:,1)]; %%% end point to the first point
% plot3(control_P(1,:),control_P(2,:),control_P(3,:),'-o'); hold on

% color = [0.6350 0.0780 0.1840];
% color = [0.4660 0.6740 0.1880];
% color = [1 0.07843 0.57647];
color = [252, 158, 127] / 255;
alf = 0.4;
% EdgeColor = [0.8549 0.43922 0.83922];
EdgeColor = [0.4940 0.1840 0.5560];
Edgealf = 0.4;
LineWidth = 0.1;
LineStyle = '-';

Bezier_polygon = Polyhedron(control_P');
B_p = plot(Bezier_polygon,'color',color); hold on

alpha(B_p,alf);
set(B_p,'LineWidth',LineWidth,'LineStyle',LineStyle,'EdgeColor',EdgeColor,'EdgeAlpha',Edgealf);
end

% function plot_piece_Bezier_polygon(A_control_P,c_control_P,zd_new,m,bar_nu_z)
% control_P = zeros(m,bar_nu_z);
% for i = 1:1:m
%     control_P(i,:) = zd_new' * A_control_P{i} + c_control_P{i};
%     % disp(control_P);
% end
% plot3(control_P(1,:),control_P(2,:),control_P(3,:),'-o'); hold on
% end


function [trajectoryfig] = plot_piece_trajectory(new_parameters,nu,t_start,t_end)
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
color = [213, 33, 32] / 255;  %%%% 红
% color = [155 34 39]/255;
% color = [0 159 105]/ 255; %%%% 绿
% color = [77 174 74] / 255; %%% 绿
% color = 'g';
trajectoryfig = plot3(Points(1,:),Points(2,:),Points(3,:),'LineWidth',0.5,'Color',color); hold on
end


function [ret] = polynomial_vector(order,t)
%%%%%% order = 多项式的最高次 + 1，0 就是常数
ret = zeros(1, order);
for i = 1:1:order
    ret(1,i) = t^(i - 1) / factorial(i - 1);
end
end

% function [ret] = reorganize_each_flatoutput(nu,W)   
% %%%% compute reorganize matrix for each flat output
% %%%% ret{i,1} * old_parameter = new_parameter
% %%%%%% nu 向量
% %%%%%% nu 向量， W 权值矩阵 W^\prime = -1/2 W，输入的 W 按照 W^\prime 使用
% W = - 0.5 * inv(W);
% 
% m = size(nu,1);
% sum_nu = sum(nu);
% 
% bar_nu_z = 2 * max(nu);
% 
% ret = cell(m,1);  %%% reorganize matrix of each flat output
% 
% raw_base_i = 0;
% for i = 1:1:m
%     reorganize_matrix_i = zeros(bar_nu_z,2*sum_nu);
%     reorganize_matrix_i(1:nu(i),raw_base_i+1:raw_base_i+nu(i)) = eye(nu(i));
%     raw_base_i = raw_base_i + nu(i);
%     col_base_j = sum_nu;
%     for j = 1:1:m
%         temp = reorganize_matrix_i(nu(i)+1:nu(i)+nu(j),col_base_j+1:col_base_j+nu(j));
%         reorganize_matrix_i(nu(i)+1:nu(i)+nu(j),col_base_j+1:col_base_j+nu(j)) = temp + W(i,j)*eye(nu(j));
%         col_base_j = col_base_j + nu(j);
%     end
%     ret{i,1} = reorganize_matrix_i;
% end
% 
% end
