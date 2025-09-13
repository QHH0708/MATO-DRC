function [tr_figs,Bp_figs] = test_b_to_trajectory_MultiAgent(nu,d,W,times,N,z0,zN,zd,tau,Reachset_c,Reachset_G,zono_or_ellipsoid)
%%% zd : a matrix (N-1) row  m*(d+1) col, interior-point constraints of trajectorys

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
bold_b_new = zd_matrix_bold_b * zd_new + vector_z0N;  %%% only applicable when d = 0
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

BP_colors = cell(N,1); %%%% Bezier polygon, Each segment has a different color
Edg_colors = cell(N,1); %%%% polygon Edgecolor;
Tr_colors = cell(N,1);   %%% trajectory color


Tr_colors{1,1} = [213, 33, 32] / 255;
Tr_colors{2,1} = [145, 60, 137] / 255;
Tr_colors{3,1} = [49, 124, 183] / 255;
Tr_colors{4,1} = [27, 124, 61] / 255;


BP_colors{1,1} = [252, 158, 127] / 255;
Edg_colors{1,1} = [0.4940 0.1840 0.5560];
BP_colors{2,1} = [246, 103, 160] / 255;
Edg_colors{2,1} = [194, 144, 197] / 255;
BP_colors{3,1} = [182, 215, 232] / 255;
Edg_colors{3,1} = [109, 173, 209] / 255;
BP_colors{4,1} = [178, 225, 185] / 255;
Edg_colors{4,1} = [119, 202, 197] / 255;


tr_figs = cell(N,1);
Bp_figs = cell(N,1);
for n = 1:1:N   %%% trajectory piece
    % bold_c_n = bold_c_all((n-1)*2*sum_nu + 1:n*2*sum_nu,1);
    bold_c_n = select_cn{n,1} * bold_c_all;
    new_parameters = zeros(m,2*max(nu));
    for i = 1:1:m  %%% for each output / each dimension
        new_parameters(i,:) = bold_c_n' * reorganize_matrices{i,1}';
    end
    % disp(new_parameters);
    r1 = (times(n) - times(1)) / tau + 1;
    r2 = (times(n+1) - times(1)) / tau;

    tr_fig = 0;
    tr_fig = plot_piece_trajectory(new_parameters,nu,times(n),times(n+1),tau,Reachset_c(r1:r2,1),Reachset_G(r1:r2,1),Tr_colors{n,1},zono_or_ellipsoid);

    %%%%%%% plot Bezier polygon
    % Bp_fig = 0;
    % Bp_fig = plot_piece_Bezier_polygon(new_parameters,bar_nu_z,times(n),times(n+1),BP_colors{n,1},Edg_colors{n,1});
    % tr_figs{n,1} = tr_fig;
    % Bp_figs{n,1} = Bp_fig;

    % trcolor = [49, 124, 183] / 255;
    % BPcolor = [182, 215, 232] / 255;
    % Edgcolor = [109, 173, 209] / 255;
    % plot_piece_trajectory(new_parameters,nu,times(n),times(n+1),tau,Reachset_c(r1:r2,1),Reachset_G(r1:r2,1),color_start,color_end,dcolor,trcolor);
    % plot_piece_Bezier_polygon(new_parameters,bar_nu_z,times(n),times(n+1),BPcolor,Edgcolor);
    
    A_CP = cell(1,m);
    c_CP = cell(1,m);
    for i = 1:1:m
        A_CP{i} = A_control_P{n,i};
        c_CP{i} = c_control_P{n,i};
    end
    % plot_piece_Bezier_polygon(A_CP,c_CP,zd_new,m,bar_nu_z);
end

end

function [B_p] = plot_piece_Bezier_polygon(new_parameters,bar_nu_z,t_start,t_end,BP_color,Edg_color)
[Ts] = Ts_time_scaling_matrix_Bernstein(bar_nu_z,t_start,t_end);
[QB] = QB_expanded_matrices_Bernstein(bar_nu_z);
control_P = new_parameters * Ts' / QB;
% disp(control_P)
% control_P = [control_P control_P(:,1)]; %%% end point to the first point
% plot3(control_P(1,:),control_P(2,:),control_P(3,:),'-o'); hold on

alf = 0.4;
Edgealf = 0.8;
LineWidth = 0.2;
LineStyle = '-';

color = BP_color;
EdgeColor = Edg_color;

control_P = control_P + 1e-6 * randn(size(control_P)); 
%%% add randm noise to the control point to address the issue of numerical
%%% accuracy when plotting Bezier polygons

Bezier_polygon = Polyhedron(control_P');
B_p = plot(Bezier_polygon,'color',color); hold on

alpha(B_p,alf);
set(B_p,'LineWidth',LineWidth,'LineStyle',LineStyle,'EdgeColor',EdgeColor,'EdgeAlpha',Edgealf);
end



function [tr_fig] = plot_piece_trajectory(new_parameters,nu,t_start,t_end,tau,Reachset_c,Reachset_G,trcolor,zono_or_ellipsoid)
%%% new_parameters is a matrix, each raw is a parameter for a dimension
%%% nus is a vector, maximal order + 1 of each dimension
%%% t_start and t_end should be the integral multiples of 0.1
dim = size(new_parameters,1);
bar_nu_z = 2 * max(nu);
num_p = (t_end - t_start) / tau + 1;
Points = zeros(dim,num_p);
point = zeros(dim,1);
accuracy = 12;
LineWidth = 0.2;
Edgealf = 0.8;
for i = 1:1:num_p
    t = t_start + (i-1) * tau;  
    for j = 1:1:dim
        point(j,1) = new_parameters(j,:) * polynomial_vector(bar_nu_z,t)';
    end
    Points(:,i) = point;

    %%% plot the reachable sets
    color = color_t(t);
    Rset_color = color;
    Rset_Edgecolor = color;
    if mod(i,10) == 1 && i <= num_p - 1
        if zono_or_ellipsoid == 0 %%% plot ellipsoid reachable sets
            plot_ellipsoid_3D(Reachset_c{i,1}+point,Reachset_G{i,1},Rset_color,0.4,LineWidth,'-',Rset_Edgecolor,Edgealf,accuracy); hold on
        end
        if zono_or_ellipsoid == 1 %%% plot zonotope reachable sets
            plot_constrained_zonotope_3D(Reachset_c{i,1}+point,Reachset_G{i,1},[],[],Rset_color,0.4,LineWidth,'-',Rset_Edgecolor,Edgealf); hold on
        end
    end
end

color = trcolor;  %%% trajectory color
tr_fig = 0;
%%%% plot the trajectory
% tr_fig = plot3(Points(1,:),Points(2,:),Points(3,:),'LineWidth',1,'Color',color); hold on

end


function [ret] = polynomial_vector(order,t)
%%%%%% order = The highest power of a polynomial + 1，
ret = zeros(1, order);
for i = 1:1:order
    ret(1,i) = t^(i - 1) / factorial(i - 1);
end
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

function [colorbar] = color_bar()
%%%% gradient color
dt = 0.1;
t0 = 0;
tN = 20;
ni = tN / dt;
colorbar = zeros(ni+1,3); % 
for i = 0:1:ni
    t = t0 + i * dt;
    colorbar(i+1,:) = color_t(t);
end
end