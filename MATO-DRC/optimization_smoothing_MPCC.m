function [Solution_num,Solutions,obj_values,exitflags,computing_time] = optimization_smoothing_MPCC(G_ellip,c_ellip,G_zono,c_zono,G_reach,c_reach,nu,d,W,times,N,z0,zN,try_times)
%%%% G_ellip,c_ellip, G_zono, c_zono are cell types
%%%% reachableSets cell{N * num_agent} reachable sets of agent,
%%%% num_agent，agent number
%%%% N：trajectory segments

% Decompose the generalized second-order cone into multiple independent ordinary second-order cones
% and zonotopes
% Arrange variables according to their disassembled form,
% as well as Lagrange multiplier variables

bar_nu_z = 2*max(nu); %%% polynomial degree + 1
m = size(nu,1); %%% flat output number
sum_nu = sum(nu);

[bold_M,zd_matrix_bold_b,vector_z0N] = bold_M_matrix(nu,d,W,times,N,z0,zN);
%%% bold_b = matrix_zd_bold_b * zd + vector_z0N , zd: variables
Wint = Wint_matrix_ofobjectivefunction(nu,W,times,N);
MtWM = bold_M' \ Wint / bold_M;

%%%% objective value is transformed to: obj = zd' * A_obj * zd + B_obj * zd + c_obj
A_obj = zd_matrix_bold_b' *  MtWM * zd_matrix_bold_b;
B_obj = 2 * vector_z0N' * MtWM * zd_matrix_bold_b;
c_obj = vector_z0N' * MtWM * vector_z0N;

%%%%%  bold_c_n = select_cn{i,1} * bold_c_all 
select_cn = cell(N,1);
for i = 1:1:N
    select_cn_temp = zeros(2*sum_nu,2*N*sum_nu);
    select_cn_temp(:,(i-1)*2*sum_nu + 1:i*2*sum_nu) = eye(2*sum_nu);
    select_cn{i,1} = select_cn_temp;
end

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
        A_control_P{n,i} = zd_matrix_bold_b' / bold_M' * select_cn{n}' * reorganize_matrices{i}' * TsQB{n};
        c_control_P{n,i} = vector_z0N' / bold_M' * select_cn{n}' * reorganize_matrices{i}' * TsQB{n};
    end
end

%%%% variable number
dim_x = (N-1)*m*(d+1) + (37*2 + 39*2 ) * N; 
dim_nlcon = ((29+9)*2 + (31+9)*2) * N; %%%% dimension of constraints (including linear variables)

%%%% variable low and upper bound: lb <= x <= ub
lb = -ones(dim_x,1)*Inf;
ub = ones(dim_x,1)*Inf;
% for i = 1:1:(N-1)*m*(d+1)
%     lb(i) = -Inf;
%     ub(i) = Inf;
% end
bx = (N-1)*m*(d+1); 
for n = 1:1:N
    for i = 1:1:2  %%% 2 Ellipsoid
        lb(bx+1) = 1;
        lb(bx+8:bx+13) = zeros(6,1);
        ub(bx+8:bx+13) = ones(6,1);
        lb(bx+14) = 0;
        lb(bx+18) = 0;
        lb(bx+22:bx+27) = zeros(6,1);
        lb(bx+28:bx+33) = zeros(6,1);
        bx = bx + 37;
    end
    for i = 1:1:2 %%%%2 zonotope
        lb(bx+1) = 1;
        lb(bx+8:bx+13) = zeros(6,1);
        ub(bx+8:bx+13) = ones(6,1);
        lb(bx+14:bx+16) = zeros(3,1);
        lb(bx+17:bx+19) = zeros(3,1);
        lb(bx+20) = 0;
        lb(bx+24:bx+29) = zeros(6,1);
        lb(bx+30:bx+35) = zeros(6,1);
        bx = bx + 39;
    end
end

%%%% varep: small threshold for scaling coefficient
varep = 0.0001;

%%%% constraint low and upper bound
lb_c = zeros(dim_nlcon,1);
ub_c= zeros(dim_nlcon,1);
bc = 0;
for n = 1:1:N
    for i = 1:1:2  %%%% 2 Ellipsoid
        lb_c(bc+1:bc+3) = - (c_ellip{i} + c_reach{n});
        ub_c(bc+1:bc+3) = - (c_ellip{i} + c_reach{n});
        lb_c(bc+30) = 1 + varep;
        ub_c(bc+30) = Inf;
        lb_c(bc+31) = 1;
        ub_c(bc+31) = 1;
        lb_c(bc+32) = 1;
        ub_c(bc+32) = 1;
        bc = bc + 38;
    end
    for i = 1:1:2  %%%% 2 zonotope
        lb_c(bc+1:bc+3) = - (c_zono{i} + c_reach{n});
        ub_c(bc+1:bc+3) = - (c_zono{i} + c_reach{n});
        lb_c(bc+32) = 1 + varep;
        ub_c(bc+32) = Inf;
        lb_c(bc+33) = 1;
        ub_c(bc+33) = 1;
        lb_c(bc+34) = 1;
        ub_c(bc+34) = 1;
        bc = bc + 40;
    end
end

%%%% generate the structures of the 1st 2nd gradient matrices
ep = 1;
epp = 1;
[jacstr_static] = jac_structure_generation(dim_x,dim_nlcon,epp,A_control_P,c_control_P,G_ellip,G_zono,G_reach,d,m,N,bar_nu_z);
[hessstr_static] = hessian_structure_generation(dim_x,dim_nlcon,A_obj,A_control_P,ep,d,m,N);


ipoptset('linear_solver', 'ma57'); %% set the linear solver for ipopt: ma57, pardiso, mumps
options = optiset('solver','ipopt','display','off'); % maxiter',1500
% options = optiset('solver','ipopt','display','iter','derivCheck','on'); %%%  derivCheck: Check gradient information


dx_num = (N-1)*m*(d+1);
Solution_num = 0;
Solutions = zeros(dx_num,try_times);
obj_values = zeros(1,try_times);
exitflags = zeros(1,try_times);
% computing_times = zeros(1,try_times);
rng(1);   %%%% random seed
%%%%%% try_times:  Smoothing MPCC will be repeatedly solved multiple times
ep_0 = 1;  %%% initial smoothing parameter

t1 = datetime('now');
for try_time = 1:1:try_times
    disp(try_time);
    %%%% random initialization of variable values
    x0 = zeros(dim_x,1);
    x0(1:(N-1)*m*(d+1)) = 50 * (rand((N-1)*m*(d+1),1) - 0.2);

    bx = (N-1)*m*(d+1);

    for n = 1:1:N  %%% N segements
        for i = 1:1:2 %%% correspond to two Ellipsoid obstacles
            x0(bx+1) = 7 * rand(1,1);
            x0(bx+2:bx+4) = 8 * (rand(3,1) - 0.5);
            x0(bx+5:bx+7) = 8 * (rand(3,1) - 0.5);
            x0(bx+8:bx+13) = rand(6,1);
            x0(bx+14) = 7 * rand(1,1);
            x0(bx+15:bx+17) = 8 * (rand(3,1) - 0.5);
            x0(bx+18) = 7 * rand(1,1);
            x0(bx+19:bx+21) = 8 * (rand(3,1) - 0.5);
            x0(bx+22:bx+27) = 8 * rand(6,1);
            x0(bx+28:bx+33) = 8 * rand(6,1);
            x0(bx+34:bx+36) = 12 * (rand(3,1) - 0.5);
            x0(bx+37) = 12 * (rand(1,1) - 0.5);
            bx = bx + 37;
        end
        for i = 1:1:2 %%% correspond to two Zonotope obstacles
            x0(bx+1) = 7 * rand(1,1);
            x0(bx+2:bx+4) = 8 * (rand(3,1) - 0.5);
            x0(bx+5:bx+7) = 8 * (rand(3,1) - 0.5);
            x0(bx+8:bx+13) = rand(6,1);
            x0(bx+14:bx+16) = 12 * rand(3,1);
            x0(bx+17:bx+19) = 12 * rand(3,1);
            x0(bx+20) = 7 * rand(1,1);
            x0(bx+21:bx+23) = 8 * (rand(3,1) - 0.5);
            x0(bx+24:bx+29) = 8 * rand(6,1);
            x0(bx+30:bx+35) = 8 * rand(6,1);
            x0(bx+36:bx+38) = 12 * (rand(3,1) - 0.5);
            x0(bx+39) = 12 * (rand(1,1) - 0.5);
            bx = bx + 39;
        end
    end  %%%% the end of random initialization of variables
   
    ep = ep_0;
    xk = x0;
    while(1)
        epp = 2*ep^2;
        %%%% Ipopt has an option to approximate the Hessian of the Lagrangian by a limited-memory quasi-Newton method (L-BFGS)
        %%%% For small problem size, LBFGS is more efficient.
        if N <= 4 %%% Do not use exact hessian matrix in small problem size.
            optimization_problem = opti('fun',@(x)objFun(x,A_obj,B_obj,c_obj,d,m,N),'grad',@(x)grad_objFun(x,dim_x,A_obj,B_obj,d,m,N),'nl',@(x)nonlinearConstraints(x,dim_nlcon,epp,A_control_P,c_control_P,G_ellip,G_zono,G_reach,d,m,N,bar_nu_z),lb_c,ub_c,'jac',@(x)jac_nlcon(x,dim_x,dim_nlcon,epp,A_control_P,c_control_P,G_ellip,G_zono,G_reach,d,m,N,bar_nu_z),'jacstr',@()jac_sparse_structure(jacstr_static),'bounds',lb,ub,'x0',xk,'options',options);
        else  %%% Use exact hessian matrix in large problem size
            optimization_problem = opti('fun',@(x)objFun(x,A_obj,B_obj,c_obj,d,m,N),'grad',@(x)grad_objFun(x,dim_x,A_obj,B_obj,d,m,N),'nl',@(x)nonlinearConstraints(x,dim_nlcon,epp,A_control_P,c_control_P,G_ellip,G_zono,G_reach,d,m,N,bar_nu_z),lb_c,ub_c,'jac',@(x)jac_nlcon(x,dim_x,dim_nlcon,epp,A_control_P,c_control_P,G_ellip,G_zono,G_reach,d,m,N,bar_nu_z),'jacstr',@()jac_sparse_structure(jacstr_static),'hess',@(x,sigma,lambda)hessian(x,sigma,lambda,dim_x,A_obj,A_control_P,ep,d,m,N),'hstr',@()hessian_structure(hessstr_static),'bounds',lb,ub,'x0',xk,'options',options);
        end
        % optimization_problem = opti('fun',@(x)objFun(x,A_obj,B_obj,c_obj,d,m,N),'grad',@(x)grad_objFun(x,dim_x,A_obj,B_obj,d,m,N),'nl',@(x)nonlinearConstraints(x,dim_nlcon,epp,A_control_P,c_control_P,G_ellip,G_zono,G_reach,d,m,N,bar_nu_z),lb_c,ub_c,'jac',@(x)jac_nlcon(x,dim_x,dim_nlcon,epp,A_control_P,c_control_P,G_ellip,G_zono,G_reach,d,m,N,bar_nu_z),'jacstr',@()jac_sparse_structure(jacstr_static),'hess',@(x,sigma,lambda)hessian(x,sigma,lambda,dim_x,A_obj,A_control_P,ep,d,m,N),'hstr',@()hessian_structure(hessstr_static),'bounds',lb,ub,'x0',xk,'options',options);
        % optimization_problem = opti('fun',@(x)objFun(x,A_obj,B_obj,c_obj,d,m,N),'grad',@(x)grad_objFun(x,dim_x,A_obj,B_obj,d,m,N),'nl',@(x)nonlinearConstraints(x,dim_nlcon,epp,A_control_P,c_control_P,G_ellip,G_zono,G_reach,d,m,N,bar_nu_z),lb_c,ub_c,'jac',@(x)jac_nlcon(x,dim_x,dim_nlcon,epp,A_control_P,c_control_P,G_ellip,G_zono,G_reach,d,m,N,bar_nu_z),'jacstr',@()jac_sparse_structure(jacstr_static),'bounds',lb,ub,'x0',xk,'options',options);
        % optimization_problem = opti('fun',@(x)objFun(x,A_obj,B_obj,c_obj,d,m,N),'nl',@(x)nonlinearConstraints(x,dim_nlcon,ep,A_control_P,c_control_P,G_ellip,G_zono,G_reach,d,m,N,bar_nu_z),lb_c,ub_c,'bounds',lb,ub,'x0',xk,'options',options);
        [xkk,fval_kk,exitflag] = solve(optimization_problem);
        if exitflag == 1 %%%% (local) optimal solution found
            if ep <= 0.01
                Solution_num = Solution_num + 1;
                Solutions(:,try_time) = xkk(1:dx_num);
                obj_values(try_time) = fval_kk;
                exitflags(try_time) = exitflag;
                break;
            else
                ep = ep / 10;
                xk = xkk;
            end
        elseif exitflag == 2  %%%% probably (local) optimal solution found
            Solution_num = Solution_num + 1;
            Solutions(:,try_time) = xkk(1:dx_num);
            obj_values(try_time) = fval_kk;
            exitflags(try_time) = exitflag;
            break;
        else
            break;
        end
    end
end
t2 = datetime('now');
computing_time = seconds(t2 - t1);

end


function [obj] = objFun(x,A_obj,B_obj,c_obj,d,m,N)  %%% objective value
obj = x(1:(N-1)*m*(d+1))' * A_obj * x(1:(N-1)*m*(d+1)) + B_obj * x(1:(N-1)*m*(d+1)) + c_obj;
end

function [grad] = grad_objFun(x,dim_x,A_obj,B_obj,d,m,N) %%% gradient
%%%% dim_x: variable number,  (N-1)*m*(d+1): variable number for the parameters of interior-point constraints
grad = zeros(1,dim_x);
grad(1,1:(N-1)*m*(d+1)) = 2 * x(1:(N-1)*m*(d+1))' * A_obj + B_obj;
end

function [nlcon] = nonlinearConstraints(x,dim_nlcon,epp,A_control_P,c_control_P,G_ellip,G_zono,G_reach,d,m,N,bar_nu_z)  %%% 非线性约束
%%%% G_ellip = cell(2,1): ellipsoid obstacles
%%%% G_zono = cell(2,1): zonotopeobstacles
%%%% G_reach = cell(N,1): reachable set of each segement
%%%% ep: smoothing parameter
%%%% Linear constraints are also included in this function!
%%%% epp = 2*ep^2;

nlcon = zeros(dim_nlcon,1);

bc = 0;  %%%% base index for constraints
bx = (N-1)*m*(d+1);  %%%% base index for variables
 
for n = 1:1:N %%% trajectroy segments
    control_P = zeros(m,bar_nu_z);
    for i = 1:1:m
        control_P(i,:) = x(1:(N-1)*m*(d+1))' * A_control_P{n,i} + c_control_P{n,i};
    end

    for j = 1:1:2 %%%% correspond 2 ellipsoid obstacles
        nlcon(bc+1:bc+3) = G_ellip{j}*x(bx+2:bx+4) + G_reach{n}*x(bx+5:bx+7) - control_P*x(bx+8:bx+13);
        nlcon(bc+4:bc+9) = x(bx+22:bx+27) - x(bx+28:bx+33) - control_P'*x(bx+34:bx+36) + ones(6,1)*x(bx+37);

        %%% w = x^2 + y^2 + 2*ep^2, w =[w1;w2]  w1_E subscript "E" dentoes Ellipsoid
        w1_E = x(bx+14:bx+17)'*x(bx+14:bx+17) + x(bx+1:bx+4)'*x(bx+1:bx+4) + epp;
        w2_E = 2 * ( x(bx+14)*x(bx+15:bx+17) + x(bx+1)*x(bx+2:bx+4) );
        %%% the vector w's spectral values and spectral vectors
        norm_w2_E = sqrt(w2_E'*w2_E);
        lamd1_E = w1_E - norm_w2_E;
        lamd2_E = w1_E + norm_w2_E;
        if norm_w2_E == 0
            % v_rand_E = rand(3,1);
            v_rand_E = ones(3,1);
            v_rand_E = v_rand_E / sqrt(v_rand_E'*v_rand_E);
            spvec1_E = 0.5*[1;-v_rand_E];
            spvec2_E = 0.5*[1;v_rand_E];
        else
            spvec1_E = 0.5*[1;-w2_E/norm_w2_E];
            spvec2_E = 0.5*[1;w2_E/norm_w2_E];
        end
        nlcon(bc+10:bc+13) = x(bx+14:bx+17) + x(bx+1:bx+4) - (sqrt(lamd1_E)*spvec1_E + sqrt(lamd2_E)*spvec2_E);
        
        %%%% w1_R subscript "R" dentoes Reachable set
        w1_R = x(bx+18:bx+21)'*x(bx+18:bx+21) + [x(bx+1);x(bx+5:bx+7)]'*[x(bx+1);x(bx+5:bx+7)] + epp;
        w2_R = 2 * ( x(bx+18)*x(bx+19:bx+21) + x(bx+1)*x(bx+5:bx+7) );
        norm_w2_R = sqrt(w2_R'*w2_R);
        lamd1_R = w1_R - norm_w2_R;
        lamd2_R = w1_R + norm_w2_R;
        if norm_w2_R == 0
            % v_rand_R = rand(3,1);
            v_rand_R = ones(3,1);
            v_rand_R = v_rand_R / sqrt(v_rand_R'*v_rand_R);
            spvec1_R = 0.5*[1;-v_rand_R];
            spvec2_R = 0.5*[1;v_rand_R];
        else
            spvec1_R = 0.5*[1;-w2_R/norm_w2_R];
            spvec2_R = 0.5*[1;w2_R/norm_w2_R];
        end
        nlcon(bc+14:bc+17) = x(bx+18:bx+21) + [x(bx+1);x(bx+5:bx+7)] - (sqrt(lamd1_R)*spvec1_R + sqrt(lamd2_R)*spvec2_R);
        %%% smoothing linear complementarity constraints
        nlcon(bc+18:bc+23) = x(bx+22:bx+27) + ones(6,1) - x(bx+8:bx+13) - sqrt(x(bx+22:bx+27).*x(bx+22:bx+27) + x(bx+8:bx+13).*x(bx+8:bx+13)-2*x(bx+8:bx+13)+ones(6,1) + ones(6,1)*epp);
        nlcon(bc+24:bc+29) = x(bx+28:bx+33) + x(bx+8:bx+13) - sqrt(x(bx+28:bx+33).*x(bx+28:bx+33) + x(bx+8:bx+13).*x(bx+8:bx+13) + ones(6,1)*epp);
        %%%% linear constraints
        nlcon(bc+30) = x(bx+1);
        nlcon(bc+31) = ones(1,6)*x(bx+8:bx+13);
        nlcon(bc+32) = x(bx+14) + x(bx+18);
        nlcon(bc+33:bc+35) = -x(bx+15:bx+17) + G_ellip{j}'*x(bx+34:bx+36);
        nlcon(bc+36:bc+38) = -x(bx+19:bx+21) + G_reach{n}'*x(bx+34:bx+36);
        bc = bc + 29 + 9;
        bx = bx + 37;
    end

    for j = 1:1:2 %%%% 2 zonotope obstacles
        nlcon(bc+1:bc+3) = G_zono{j}*x(bx+2:bx+4) + G_reach{n}*x(bx+5:bx+7) - control_P*x(bx+8:bx+13);
        nlcon(bc+4:bc+9) = x(bx+24:bx+29) - x(bx+30:bx+35) - control_P'*x(bx+36:bx+38) + ones(6,1)*x(bx+39);
        %%% smoothing linear complementarity constraints
        nlcon(bc+10:bc+12) = x(bx+14:bx+16) + x(bx+1)-x(bx+2:bx+4) - sqrt(x(bx+14:bx+16).^2 + (x(bx+1)-x(bx+2:bx+4)).^2 + epp);
        nlcon(bc+13:bc+15) = x(bx+17:bx+19) + x(bx+1)+x(bx+2:bx+4) - sqrt(x(bx+17:bx+19).^2 + (x(bx+1)+x(bx+2:bx+4)).^2 + epp);
        %%% smoothing second-order cone complementarity constraints
        w1 = x(bx+20:bx+23)'*x(bx+20:bx+23) + x(bx+1)^2+x(bx+5:bx+7)'*x(bx+5:bx+7) + epp;
        w2 = 2 * (x(bx+20)*x(bx+21:bx+23) + x(bx+1)*x(bx+5:bx+7));
        norm_w2 = sqrt(w2'*w2);
        lamd1 = w1 - norm_w2;
        lamd2 = w1 + norm_w2;
        if norm_w2 == 0
            % v_rand = rand(3,1);
            v_rand = ones(3,1);
            v_rand = v_rand / sqrt(v_rand'*v_rand);
            spvec1 = 0.5*[1;-v_rand];
            spvec2 = 0.5*[1;v_rand];
        else
            spvec1 = 0.5*[1;-w2/norm_w2];
            spvec2 = 0.5*[1;w2/norm_w2];
        end
        nlcon(bc+16:bc+19) = x(bx+20:bx+23) + [x(bx+1);x(bx+5:bx+7)] - (sqrt(lamd1)*spvec1 + sqrt(lamd2)*spvec2);
        %%% smoothing linear complementarity constraints
        nlcon(bc+20:bc+25) = x(bx+24:bx+29) + 1-x(bx+8:bx+13) - sqrt(x(bx+24:bx+29).^2 + (1-x(bx+8:bx+13)).^2 + epp);
        nlcon(bc+26:bc+31) = x(bx+30:bx+35) + x(bx+8:bx+13) - sqrt(x(bx+30:bx+35).^2 + x(bx+8:bx+13).^2 + epp);
        %%%% linear constraints
        nlcon(bc+32) = x(bx+1);
        nlcon(bc+33) = ones(1,6)*x(bx+8:bx+13);
        nlcon(bc+34) = ones(1,3)*x(bx+14:bx+16) + ones(1,3)*x(bx+17:bx+19) + x(bx+20);
        nlcon(bc+35:bc+37) = x(bx+14:bx+16) - x(bx+17:bx+19) + G_zono{j}'*x(bx+36:bx+38);
        nlcon(bc+38:bc+40) = -x(bx+21:bx+23) + G_reach{n}'*x(bx+36:bx+38);
        bc = bc + 31 + 9;
        bx = bx + 39;
    end

end

end

function [jac] = jac_nlcon(x,dim_x,dim_nlcon,epp,A_control_P,c_control_P,G_ellip,G_zono,G_reach,d,m,N,bar_nu_z) %%% 非线性约束的梯度
%%%%% epp = 2*ep^2;

%%% Jacobian matrix of constraints, row number,
%%% row number=constraint number, col number = variable number
jac = zeros(dim_nlcon,dim_x);

bc = 0;  %%%% base index for constraints
dx = (N-1)*m*(d+1); 
bx = (N-1)*m*(d+1);  %%%% base index for variables

for n = 1:1:N %%%% N segements
    control_P = zeros(m,bar_nu_z);
    for i = 1:1:m
        control_P(i,:) = x(1:dx)' * A_control_P{n,i} + c_control_P{n,i};
    end

    for j = 1:1:2  %%%% 2 ellipsoid obstacles
        jac(bc+1:bc+3,bx+2:bx+4) = G_ellip{j};
        jac(bc+1:bc+3,bx+5:bx+7) = G_reach{n};
        jac(bc+1:bc+3,bx+8:bx+13) = -control_P;
        for i = 1:1:m
            jac(bc+i,1:dx) = - (A_control_P{n,i} * x(bx+8:bx+13))';
        end
        jac(bc+4:bc+9,bx+22:bx+27) = eye(6);
        jac(bc+4:bc+9,bx+28:bx+33) = -eye(6);
        jac(bc+4:bc+9,bx+34:bx+36) = -control_P';
        jac(bc+4:bc+9,1:dx) = - x(bx+34)*A_control_P{n,1}' - x(bx+35)*A_control_P{n,2}' - x(bx+36)*A_control_P{n,3}';
        jac(bc+4:bc+9,bx+37) = ones(6,1);
        %%%% Jacobian matrix of second-order cone complementarity
        %%%% constraints
        %%%% w = x^2 + y^2 + 2*ep^2, w =[w1;w2], subscript "E" for ellipsoid
        w1_E = x(bx+14:bx+17)'*x(bx+14:bx+17) + x(bx+1:bx+4)'*x(bx+1:bx+4) + epp;
        w2_E = 2 * ( x(bx+14)*x(bx+15:bx+17) + x(bx+1)*x(bx+2:bx+4) );
        %%% the vector w's spectral values and spectral vectors
        norm_w2_E = sqrt(w2_E'*w2_E);
        lamd1_E = w1_E - norm_w2_E;
        lamd2_E = w1_E + norm_w2_E;
        if norm_w2_E == 0
            % v_rand_E = rand(3,1);
            v_rand_E = ones(3,1);
            v_rand_E = v_rand_E / sqrt(v_rand_E'*v_rand_E);
            spvec1_E = 0.5*[1;-v_rand_E];
            spvec2_E = 0.5*[1;v_rand_E];
        else
            spvec1_E = 0.5*[1;-w2_E/norm_w2_E];
            spvec2_E = 0.5*[1;w2_E/norm_w2_E];
        end
        sqrtw_E = sqrt(lamd1_E)*spvec1_E + sqrt(lamd2_E)*spvec2_E;
        L_sqrtw_E = [sqrtw_E';sqrtw_E(2:4) sqrtw_E(1)*eye(3)];
        L_lamda_E = [x(bx+14:bx+17)';x(bx+15:bx+17) x(bx+14)*eye(3)];
        L_rxi_E = [x(bx+1:bx+4)';x(bx+2:bx+4) x(bx+1)*eye(3)];
        jac(bc+10:bc+13,bx+14:bx+17) = (eye(4) - L_lamda_E / L_sqrtw_E)';
        jac(bc+10:bc+13,bx+1:bx+4) = (eye(4) - L_rxi_E / L_sqrtw_E)';
        %%%% w1_R subscript "R" for Reachable set
        w1_R = x(bx+18:bx+21)'*x(bx+18:bx+21) + [x(bx+1);x(bx+5:bx+7)]'*[x(bx+1);x(bx+5:bx+7)] + epp;
        w2_R = 2 * ( x(bx+18)*x(bx+19:bx+21) + x(bx+1)*x(bx+5:bx+7) );
        norm_w2_R = sqrt(w2_R'*w2_R);
        lamd1_R = w1_R - norm_w2_R;
        lamd2_R = w1_R + norm_w2_R;
        if norm_w2_R == 0
            % v_rand_R = rand(3,1);
            v_rand_R = ones(3,1);
            v_rand_R = v_rand_R / sqrt(v_rand_R'*v_rand_R);
            spvec1_R = 0.5*[1;-v_rand_R];
            spvec2_R = 0.5*[1;v_rand_R];
        else
            spvec1_R = 0.5*[1;-w2_R/norm_w2_R];
            spvec2_R = 0.5*[1;w2_R/norm_w2_R];
        end
        sqrtw_R = sqrt(lamd1_R)*spvec1_R + sqrt(lamd2_R)*spvec2_R;
        L_sqrtw_R = [sqrtw_R';sqrtw_R(2:4) sqrtw_R(1)*eye(3)];
        L_lamda_R = [x(bx+18:bx+21)';x(bx+19:bx+21) x(bx+18)*eye(3)];
        L_rxi_R = [x(bx+1) x(bx+5:bx+7)';x(bx+5:bx+7) x(bx+1)*eye(3)];
        jac(bc+14:bc+17,bx+18:bx+21) = (eye(4) - L_lamda_R / L_sqrtw_R)';
        temp = (eye(4) - L_rxi_R / L_sqrtw_R)';
        jac(bc+14:bc+17,bx+1) = temp(:,1);
        jac(bc+14:bc+17,bx+5:bx+7) = temp(:,2:4);
        %%%% Jacobian for linear complementarity constraints
        atemp1 = sqrt(x(bx+22:bx+27).^2 + (1-x(bx+8:bx+13)).^2 + epp);
        jac(bc+18:bc+23,bx+22:bx+27) = eye(6) - diag(x(bx+22:bx+27)./atemp1);
        jac(bc+18:bc+23,bx+8:bx+13) = -eye(6) + diag((1-x(bx+8:bx+13))./atemp1);
        atemp2 = sqrt(x(bx+28:bx+33).^2 + x(bx+8:bx+13).^2 + epp);
        jac(bc+24:bc+29,bx+28:bx+33) = eye(6) - diag(x(bx+28:bx+33)./atemp2);
        jac(bc+24:bc+29,bx+8:bx+13) = eye(6) - diag(x(bx+8:bx+13)./atemp2);
        %%%% Jacobian for linear constraints
        jac(bc+30,bx+1) = 1;
        jac(bc+31,bx+8:bx+13) = ones(1,6);
        jac(bc+32,bx+14) = 1;
        jac(bc+32,bx+18) = 1;
        jac(bc+33:bc+35,bx+15:bx+17) = - eye(3);
        jac(bc+33:bc+35,bx+34:bx+36) = G_ellip{j}';
        jac(bc+36:bc+38,bx+19:bx+21) = - eye(3);
        jac(bc+36:bc+38,bx+34:bx+36) = G_reach{n}';
        bc = bc + 38;
        bx = bx + 37;
    end

    for j = 1:1:2  %%%% 2 zonotope obstacles
        jac(bc+1:bc+3,bx+2:bx+4) = G_zono{j};
        jac(bc+1:bc+3,bx+5:bx+7) = G_reach{n};
        jac(bc+1:bc+3,bx+8:bx+13) = - control_P;
        for i = 1:1:m
            jac(bc+i,1:dx) = - (A_control_P{n,i} * x(bx+8:bx+13))';
        end
        jac(bc+4:bc+9,bx+24:bx+29) = eye(6);
        jac(bc+4:bc+9,bx+30:bx+35) = - eye(6);
        jac(bc+4:bc+9,bx+36:bx+38) = - control_P';
        jac(bc+4:bc+9,1:dx) = - x(bx+36)*A_control_P{n,1}' - x(bx+37)*A_control_P{n,2}' - x(bx+38)*A_control_P{n,3}';
        jac(bc+4:bc+9,bx+39) = ones(6,1);
        %%%% jacobian for linear complementarity constraints
        atemp1 = sqrt(x(bx+14:bx+16).^2 + (x(bx+1)-x(bx+2:bx+4)).^2 + epp);
        jac(bc+10:bc+12,bx+14:bx+16) = eye(3) - diag(x(bx+14:bx+16)./atemp1);
        jac(bc+10:bc+12,bx+1) = ones(3,1) - (x(bx+1)-x(bx+2:bx+4))./atemp1;
        jac(bc+10:bc+12,bx+2:bx+4) = - eye(3) + diag((x(bx+1)-x(bx+2:bx+4))./atemp1);
        atemp2 = sqrt(x(bx+17:bx+19).^2 + (x(bx+1)+x(bx+2:bx+4)).^2 + epp);
        jac(bc+13:bc+15,bx+17:bx+19) = eye(3) - diag(x(bx+17:bx+19)./atemp2);
        jac(bc+13:bc+15,bx+1) = ones(3,1) - (x(bx+1)+x(bx+2:bx+4))./atemp2;
        jac(bc+13:bc+15,bx+2:bx+4) = eye(3) - diag((x(bx+1)+x(bx+2:bx+4))./atemp2);
        %%%% jacobian for linear complementarity constraints: reachable set
        w1 = x(bx+20:bx+23)'*x(bx+20:bx+23) + x(bx+1)^2+x(bx+5:bx+7)'*x(bx+5:bx+7) + epp;
        w2 = 2 * (x(bx+20)*x(bx+21:bx+23) + x(bx+1)*x(bx+5:bx+7));
        norm_w2 = sqrt(w2'*w2);
        lamd1 = w1 - norm_w2;
        lamd2 = w1 + norm_w2;
        if norm_w2 == 0
            % v_rand = rand(3,1);
            v_rand = ones(3,1);
            v_rand = v_rand / sqrt(v_rand'*v_rand);
            spvec1 = 0.5*[1;-v_rand];
            spvec2 = 0.5*[1;v_rand];
        else
            spvec1 = 0.5*[1;-w2/norm_w2];
            spvec2 = 0.5*[1;w2/norm_w2];
        end
        sqrtw = sqrt(lamd1)*spvec1 + sqrt(lamd2)*spvec2;
        L_sqrtw = [sqrtw';sqrtw(2:4) sqrtw(1)*eye(3)];
        L_lamdaR = [x(bx+20:bx+23)';x(bx+21:bx+23) x(bx+20)*eye(3)];
        L_rxiR = [x(bx+1) x(bx+5:bx+7)';x(bx+5:bx+7) x(bx+1)*eye(3)];
        jac(bc+16:bc+19,bx+20:bx+23) = (eye(4) - L_lamdaR / L_sqrtw)';
        temp = (eye(4) - L_rxiR / L_sqrtw)';
        jac(bc+16:bc+19,bx+1) = temp(:,1);
        jac(bc+16:bc+19,bx+5:bx+7) = temp(:,2:4);
        %%%% jacobian for linear complementarity constraints: Bezier polygon
        atemp1 = sqrt(x(bx+24:bx+29).^2 + (1-x(bx+8:bx+13)).^2 + epp);
        jac(bc+20:bc+25,bx+24:bx+29) = eye(6) - diag(x(bx+24:bx+29)./atemp1);
        jac(bc+20:bc+25,bx+8:bx+13) = -eye(6) + diag((1-x(bx+8:bx+13))./atemp1);
        atemp2 = sqrt(x(bx+30:bx+35).^2 + x(bx+8:bx+13).^2 + epp);
        jac(bc+26:bc+31,bx+30:bx+35) = eye(6) - diag(x(bx+30:bx+35)./atemp2);
        jac(bc+26:bc+31,bx+8:bx+13) = eye(6) - diag(x(bx+8:bx+13)./atemp2);
        %%%% jacobian for linear constraints
        jac(bc+32,bx+1) = 1;
        jac(bc+33,bx+8:bx+13) = ones(1,6);
        jac(bc+34,bx+14:bx+16) = ones(1,3);
        jac(bc+34,bx+17:bx+19) = ones(1,3);
        jac(bc+34,bx+20) = 1;
        jac(bc+35:bc+37,bx+14:bx+16) = eye(3);
        jac(bc+35:bc+37,bx+17:bx+19) = - eye(3);
        jac(bc+35:bc+37,bx+36:bx+38) = G_zono{j}';
        jac(bc+38:bc+40,bx+21:bx+23) = - eye(3);
        jac(bc+38:bc+40,bx+36:bx+38) = G_reach{n}';
        bc = bc + 40;
        bx = bx + 39;
    end
end
end

function [jacstr] = jac_sparse_structure(jacstr_static)
jacstr = jacstr_static;
end

function [jacstr] = jac_structure_generation(dim_x,dim_nlcon,epp,A_control_P,c_control_P,G_ellip,G_zono,G_reach,d,m,N,bar_nu_z)
x = 5 * rand(dim_x,1);
[jac] = jac_nlcon(x,dim_x,dim_nlcon,epp,A_control_P,c_control_P,G_ellip,G_zono,G_reach,d,m,N,bar_nu_z);
jac(jac ~= 0) = 1;
jacstr = sparse(jac);
end

function [hess] = hessian(x,sigma,lambda,dim_x,A_obj,A_control_P,ep,d,m,N)
dx_num = (N-1)*m*(d+1);
hess = zeros(dim_x,dim_x);  %%%%% initialization hessian matrix
hess(1:dx_num,1:dx_num) = sigma * 2 * tril(A_obj);  %%%% 2nd gradient of obj
J1 = 2*eye(4);
J2 = 2*[0 1 0 0;1 0 0 0;0 0 0 0;0 0 0 0];
J3 = 2*[0 0 1 0;0 0 0 0;1 0 0 0;0 0 0 0];
J4 = 2*[0 0 0 1;0 0 0 0;0 0 0 0;1 0 0 0];
bc = 0;
bx = dx_num;  %%% base index
for n = 1:1:N
    for j = 1:1:2  %%% 2 ellipsoid obstacle
        for i = 1:1:m  %%%% m = 3
            % hess(1:dx_num,bx+8:bx+13) = hess(1:dx_num,bx+8:bx+13) - lambda(bc+i) * A_control_P{n,i};
            hess(bx+8:bx+13,1:dx_num) = hess(bx+8:bx+13,1:dx_num) - lambda(bc+i) * A_control_P{n,i}';
        end
        for i = 1:1:m
            % hess(1:dx_num,bx+33+i) = - A_control_P{n,i} * lambda(bc+4:bc+9);
            hess(bx+33+i,1:dx_num) = - (A_control_P{n,i} * lambda(bc+4:bc+9))';
        end
        %%% w = x^2 + y^2 + 2*ep^2, w =[w1;w2]  w1_E 下标E 代表Ellipsoid
        w1_E = x(bx+14:bx+17)'*x(bx+14:bx+17) + x(bx+1:bx+4)'*x(bx+1:bx+4) + 2*ep^2;
        w2_E = 2 * ( x(bx+14)*x(bx+15:bx+17) + x(bx+1)*x(bx+2:bx+4) );
        %%% the vector w's spectral values and spectral vectors
        norm_w2_E = sqrt(w2_E'*w2_E);
        lamd1_E = w1_E - norm_w2_E;
        lamd2_E = w1_E + norm_w2_E;
        if norm_w2_E == 0
            % v_rand_E = rand(3,1);
            v_rand_E = ones(3,1);
            v_rand_E = v_rand_E / sqrt(v_rand_E'*v_rand_E);
            spvec1_E = 0.5*[1;-v_rand_E];
            spvec2_E = 0.5*[1;v_rand_E];
        else
            spvec1_E = 0.5*[1;-w2_E/norm_w2_E];
            spvec2_E = 0.5*[1;w2_E/norm_w2_E];
        end
        w = sqrt(lamd1_E)*spvec1_E + sqrt(lamd2_E)*spvec2_E;
        Lw = [w';w(2:4) w(1)*eye(3)];
        iLw = inv(Lw);
        %%%%% left for x, right for y
        Lx = [x(bx+14:bx+17)';x(bx+15:bx+17) x(bx+14)*eye(3)];
        Ly = [x(bx+1:bx+4)';x(bx+2:bx+4) x(bx+1)*eye(3)];
        temp = zeros(4,4);
        for i = 1:1:4
            temp = temp + lambda(bc+9+i) * (1/8) * ( iLw(i,1)*(2*Lx*iLw*J1*iLw*2*Lx) + iLw(i,2)*(2*Lx*iLw*J2*iLw*2*Lx) + iLw(i,3)*(2*Lx*iLw*J3*iLw*2*Lx) + iLw(i,4)*(2*Lx*iLw*J4*iLw*2*Lx));
            temp = temp - lambda(bc+9+i) * (1/2) * ( iLw(i,1)*J1 + iLw(i,2)*J2 + iLw(i,3)*J3 + iLw(i,4)*J4);
        end
        hess(bx+14:bx+17,bx+14:bx+17) = tril(temp);
        temp = zeros(4,4);
        for i = 1:1:4
            temp = temp + lambda(bc+9+i) * (1/8) * ( iLw(i,1)*(2*Ly*iLw*J1*iLw*2*Ly) + iLw(i,2)*(2*Ly*iLw*J2*iLw*2*Ly) + iLw(i,3)*(2*Ly*iLw*J3*iLw*2*Ly) + iLw(i,4)*(2*Ly*iLw*J4*iLw*2*Ly));
            temp = temp - lambda(bc+9+i) * (1/2) * ( iLw(i,1)*J1 + iLw(i,2)*J2 + iLw(i,3)*J3 + iLw(i,4)*J4);
        end
        hess(bx+1:bx+4,bx+1:bx+4) = tril(temp);
        for i = 1:1:4
            % hess(bx+1:bx+4,bx+14:bx+17) = hess(bx+1:bx+4,bx+14:bx+17) + lambda(bc+9+i) * (1/8) * ( iLw(i,1)*(2*Ly*iLw*J1*iLw*2*Lx) + iLw(i,2)*(2*Ly*iLw*J2*iLw*2*Lx) + iLw(i,3)*(2*Ly*iLw*J3*iLw*2*Lx) + iLw(i,4)*(2*Ly*iLw*J4*iLw*2*Lx) );
            hess(bx+14:bx+17,bx+1:bx+4) = hess(bx+14:bx+17,bx+1:bx+4) + lambda(bc+9+i) * (1/8) * ( iLw(i,1)*(2*Lx*iLw*J1*iLw*2*Ly) + iLw(i,2)*(2*Lx*iLw*J2*iLw*2*Ly) + iLw(i,3)*(2*Lx*iLw*J3*iLw*2*Ly) + iLw(i,4)*(2*Lx*iLw*J4*iLw*2*Ly) );
        end
        % hess(bx+14:bx+17,bx+1:bx+4) = hess(bx+1:bx+4,bx+14:bx+17)';

        %%%% w1_R, "R": Reachable set
        w1_R = x(bx+18:bx+21)'*x(bx+18:bx+21) + [x(bx+1);x(bx+5:bx+7)]'*[x(bx+1);x(bx+5:bx+7)] + 2*ep^2;
        w2_R = 2 * ( x(bx+18)*x(bx+19:bx+21) + x(bx+1)*x(bx+5:bx+7) );
        norm_w2_R = sqrt(w2_R'*w2_R);
        lamd1_R = w1_R - norm_w2_R;
        lamd2_R = w1_R + norm_w2_R;
        if norm_w2_R == 0
            % v_rand_R = rand(3,1);
            v_rand_R = ones(3,1);
            v_rand_R = v_rand_R / sqrt(v_rand_R'*v_rand_R);
            spvec1_R = 0.5*[1;-v_rand_R];
            spvec2_R = 0.5*[1;v_rand_R];
        else
            spvec1_R = 0.5*[1;-w2_R/norm_w2_R];
            spvec2_R = 0.5*[1;w2_R/norm_w2_R];
        end
        w = sqrt(lamd1_R)*spvec1_R + sqrt(lamd2_R)*spvec2_R;
        Lw = [w';w(2:4) w(1)*eye(3)];
        iLw = inv(Lw);
        %%%%% left for x, right for y
        Lx = [x(bx+18:bx+21)';x(bx+19:bx+21) x(bx+18)*eye(3)];
        Ly = [x(bx+1) x(bx+5:bx+7)';x(bx+5:bx+7) x(bx+1)*eye(3)];
        temp = zeros(4,4);
        for i = 1:1:4
            temp = temp + lambda(bc+13+i) * (1/8) * ( iLw(i,1)*(2*Lx*iLw*J1*iLw*2*Lx) + iLw(i,2)*(2*Lx*iLw*J2*iLw*2*Lx) + iLw(i,3)*(2*Lx*iLw*J3*iLw*2*Lx) + iLw(i,4)*(2*Lx*iLw*J4*iLw*2*Lx));
            temp = temp - lambda(bc+13+i) * (1/2) * ( iLw(i,1)*J1 + iLw(i,2)*J2 + iLw(i,3)*J3 + iLw(i,4)*J4);
        end
        hess(bx+18:bx+21,bx+18:bx+21) = tril(temp);
        temp = zeros(4,4);
        for i = 1:1:4
            temp = temp + lambda(bc+13+i) * (1/8) * ( iLw(i,1)*(2*Ly*iLw*J1*iLw*2*Ly) + iLw(i,2)*(2*Ly*iLw*J2*iLw*2*Ly) + iLw(i,3)*(2*Ly*iLw*J3*iLw*2*Ly) + iLw(i,4)*(2*Ly*iLw*J4*iLw*2*Ly) );
            temp = temp - lambda(bc+13+i) * (1/2) * ( iLw(i,1)*J1 + iLw(i,2)*J2 + iLw(i,3)*J3 + iLw(i,4)*J4);
        end
        hess(bx+1,bx+1) = hess(bx+1,bx+1) + temp(1,1);
        % hess(bx+1,bx+5:bx+7) = temp(1,2:4);
        hess(bx+5:bx+7,bx+1) = temp(2:4,1);
        hess(bx+5:bx+7,bx+5:bx+7) = tril(temp(2:4,2:4));
        temp = zeros(4,4);
        for i = 1:1:4
            temp = temp + lambda(bc+13+i) * (1/8) * ( iLw(i,1)*(2*Ly*iLw*J1*iLw*2*Lx) + iLw(i,2)*(2*Ly*iLw*J2*iLw*2*Lx) + iLw(i,3)*(2*Ly*iLw*J3*iLw*2*Lx) + iLw(i,4)*(2*Ly*iLw*J4*iLw*2*Lx) );
        end
        % hess(bx+1,bx+18:bx+21) = temp(1,1:4);
        hess(bx+18:bx+21,bx+1) = temp(1,1:4)';
        % hess(bx+5:bx+7,bx+18:bx+21) = temp(2:4,1:4);
        hess(bx+18:bx+21,bx+5:bx+7) = temp(2:4,1:4)';
        %%%%% linear complementarity constraints
        atemp = (x(bx+22:bx+27).^2 + (1 - x(bx+8:bx+13)).^2 + 2*ep^2).^1.5;
        hess(bx+22:bx+27,bx+22:bx+27) = diag( lambda(bc+18:bc+23).*  ( - ((1-x(bx+8:bx+13)).^2 + 2*ep^2)./ atemp ) );
        hess(bx+8:bx+13,bx+8:bx+13) = diag( lambda(bc+18:bc+23).* ( - (x(bx+22:bx+27).^2 + 2*ep^2) ./ atemp ) );
        hess(bx+22:bx+27,bx+8:bx+13) = diag( lambda(bc+18:bc+23).* ( - x(bx+22:bx+27).*(1-x(bx+8:bx+13)) ./ atemp ) );
        % hess(bx+22:bx+27,bx+8:bx+13) = hess(bx+8:bx+13,bx+22:bx+27)';

        atemp = (x(bx+28:bx+33).^2 + x(bx+8:bx+13).^2 + 2*ep^2).^1.5;
        hess(bx+28:bx+33,bx+28:bx+33) = diag( lambda(bc+24:bc+29).* (- (x(bx+8:bx+13).^2 + 2*ep^2)./ atemp ) );
        hess(bx+8:bx+13,bx+8:bx+13) = hess(bx+8:bx+13,bx+8:bx+13) + diag( lambda(bc+24:bc+29).* (- (x(bx+28:bx+33).^2 + 2*ep^2)./ atemp));
        hess(bx+28:bx+33,bx+8:bx+13) = diag( lambda(bc+24:bc+29).* ( x(bx+8:bx+13).*x(bx+28:bx+33) ./ atemp ) );
        % hess(bx+28:bx+33,bx+8:bx+13) = hess(bx+8:bx+13,bx+28:bx+33)';

        bx = bx + 37;
        bc = bc + 38;
    end

    for j = 1:1:2 %%% 2 zonotope obstacles
        for i = 1:1:m  %%%% m = 3
            % hess(1:dx_num,bx+8:bx+13) = hess(1:dx_num,bx+8:bx+13) - lambda(bc+i) * A_control_P{n,i};
            hess(bx+8:bx+13,1:dx_num) = hess(bx+8:bx+13,1:dx_num) - lambda(bc+i) * A_control_P{n,i}';
        end
        for i = 1:1:m
            % hess(1:dx_num,bx+35+i) = - A_control_P{n,i} * lambda(bc+4:bc+9);
            hess(bx+35+i,1:dx_num) = - (A_control_P{n,i} * lambda(bc+4:bc+9))';
        end
        %%%% linear complementarity constraints
        atemp = ( x(bx+14:bx+16).^2 + (x(bx+1) - x(bx+2:bx+4)).^2 + 2*ep^2 ).^1.5;
        hess(bx+14:bx+16,bx+14:bx+16) = diag( lambda(bc+10:bc+12).* (- ((x(bx+1) - x(bx+2:bx+4)).^2 + 2*ep^2)./ atemp ) );
        hess(bx+2:bx+4,bx+2:bx+4) = diag( lambda(bc+10:bc+12).* (- (x(bx+14:bx+16).^2 +2*ep^2)./ atemp ) );
        % hess(bx+2:bx+4,bx+14:bx+16) = diag( lambda(bc+10:bc+12).* (- (x(bx+1) - x(bx+2:bx+4)).*x(bx+14:bx+16)./ atemp ) );
        hess(bx+14:bx+16,bx+2:bx+4) = diag( lambda(bc+10:bc+12).* (- (x(bx+1) - x(bx+2:bx+4)).*x(bx+14:bx+16)./ atemp ) );
        hess(bx+1,bx+1) = hess(bx+1,bx+1) + lambda(bc+10:bc+12)' * (- (x(bx+14:bx+16).^2 +2*ep^2)./ atemp);
        hess(bx+14:bx+16,bx+1) = lambda(bc+10:bc+12) .* ((x(bx+1) - x(bx+2:bx+4)).*x(bx+14:bx+16)./ atemp);
        % hess(bx+1,bx+14:bx+16) = hess(bx+14:bx+16,bx+1)';
        hess(bx+2:bx+4,bx+1) = lambda(bc+10:bc+12) .* ( (x(bx+14:bx+16).^2 + 2*ep^2)./atemp );
        % hess(bx+1,bx+2:bx+4) = hess(bx+2:bx+4,bx+1)';

        atemp = ( x(bx+17:bx+19).^2 + (x(bx+1) + x(bx+2:bx+4)).^2 + 2*ep^2 ).^1.5;
        hess(bx+17:bx+19,bx+17:bx+19) = diag( lambda(bc+13:bc+15).* (- ((x(bx+1) + x(bx+2:bx+4)).^2 + 2*ep^2)./atemp ) );
        hess(bx+2:bx+4,bx+2:bx+4) = hess(bx+2:bx+4,bx+2:bx+4) + diag( lambda(bc+13:bc+15).* (- (x(bx+17:bx+19).^2 + 2*ep^2)./atemp ) );
        % hess(bx+2:bx+4,bx+17:bx+19) = diag(lambda(bc+13:bc+15).* (x(bx+17:bx+19).* (x(bx+1) + x(bx+2:bx+4))./ atemp ) );
        hess(bx+17:bx+19,bx+2:bx+4) = diag(lambda(bc+13:bc+15).* (x(bx+17:bx+19).* (x(bx+1) + x(bx+2:bx+4))./ atemp ) );
        hess(bx+1,bx+1) = hess(bx+1,bx+1) + lambda(bc+13:bc+15)' * ( - (x(bx+17:bx+19).^2 + 2*ep^2)./ atemp );
        hess(bx+17:bx+19,bx+1) = lambda(bc+13:bc+15) .* ((x(bx+1) + x(bx+2:bx+4)).*x(bx+17:bx+19)./ atemp);
        % hess(bx+1,bx+17:bx+19) = hess(bx+17:bx+19,bx+1)';
        hess(bx+2:bx+4,bx+1) = hess(bx+2:bx+4,bx+1) + lambda(bc+13:bc+15) .* (- (x(bx+17:bx+19).^2 + 2*ep^2)./ atemp );
        % hess(bx+1,bx+2:bx+4) = hess(bx+2:bx+4,bx+1)';

        %%%%% second-order cone complementarity constraints
        w1 = x(bx+20:bx+23)'*x(bx+20:bx+23) + x(bx+1)^2+x(bx+5:bx+7)'*x(bx+5:bx+7) + 2*ep^2;
        w2 = 2 * (x(bx+20)*x(bx+21:bx+23) + x(bx+1)*x(bx+5:bx+7));
        norm_w2 = sqrt(w2'*w2);
        lamd1 = w1 - norm_w2;
        lamd2 = w1 + norm_w2;
        if norm_w2 == 0
            % v_rand = rand(3,1);
            v_rand = ones(3,1);
            v_rand = v_rand / sqrt(v_rand'*v_rand);
            spvec1 = 0.5*[1;-v_rand];
            spvec2 = 0.5*[1;v_rand];
        else
            spvec1 = 0.5*[1;-w2/norm_w2];
            spvec2 = 0.5*[1;w2/norm_w2];
        end
        w = sqrt(lamd1)*spvec1 + sqrt(lamd2)*spvec2;
        Lw = [w';w(2:4) w(1)*eye(3)];
        iLw = inv(Lw);
        %%%%% left for x, right for y
        Lx = [x(bx+20:bx+23)';x(bx+21:bx+23) x(bx+20)*eye(3)];
        Ly = [x(bx+1) x(bx+5:bx+7)';x(bx+5:bx+7) x(bx+1)*eye(3)];
        temp = zeros(4,4);
        for i = 1:1:4
            temp = temp + lambda(bc+15+i) * (1/8) * ( iLw(i,1)*(2*Lx*iLw*J1*iLw*2*Lx) + iLw(i,2)*(2*Lx*iLw*J2*iLw*2*Lx) + iLw(i,3)*(2*Lx*iLw*J3*iLw*2*Lx) + iLw(i,4)*(2*Lx*iLw*J4*iLw*2*Lx));
            temp = temp - lambda(bc+15+i) * (1/2) * ( iLw(i,1)*J1 + iLw(i,2)*J2 + iLw(i,3)*J3 + iLw(i,4)*J4);
        end
        hess(bx+20:bx+23,bx+20:bx+23) = tril(temp);
        temp = zeros(4,4);
        for i = 1:1:4
            temp = temp + lambda(bc+15+i) * (1/8) * ( iLw(i,1)*(2*Ly*iLw*J1*iLw*2*Ly) + iLw(i,2)*(2*Ly*iLw*J2*iLw*2*Ly) + iLw(i,3)*(2*Ly*iLw*J3*iLw*2*Ly) + iLw(i,4)*(2*Ly*iLw*J4*iLw*2*Ly) );
            temp = temp - lambda(bc+15+i) * (1/2) * ( iLw(i,1)*J1 + iLw(i,2)*J2 + iLw(i,3)*J3 + iLw(i,4)*J4);
        end
        hess(bx+1,bx+1) = hess(bx+1,bx+1) + temp(1,1);
        % hess(bx+1,bx+5:bx+7) = temp(1,2:4);
        hess(bx+5:bx+7,bx+1) = temp(2:4,1);
        hess(bx+5:bx+7,bx+5:bx+7) = tril(temp(2:4,2:4));
        temp = zeros(4,4);
        for i = 1:1:4
            temp = temp + lambda(bc+15+i) * (1/8) * ( iLw(i,1)*(2*Ly*iLw*J1*iLw*2*Lx) + iLw(i,2)*(2*Ly*iLw*J2*iLw*2*Lx) + iLw(i,3)*(2*Ly*iLw*J3*iLw*2*Lx) + iLw(i,4)*(2*Ly*iLw*J4*iLw*2*Lx) );
        end
        % hess(bx+1,bx+20:bx+23) = temp(1,1:4);
        hess(bx+20:bx+23,bx+1) = temp(1,1:4)';
        % hess(bx+5:bx+7,bx+20:bx+23) = temp(2:4,1:4);
        hess(bx+20:bx+23,bx+5:bx+7) = temp(2:4,1:4)';
        
        %%%%% 
        atemp = (x(bx+24:bx+29).^2 + (1 - x(bx+8:bx+13)).^2 + 2*ep^2).^1.5;
        hess(bx+24:bx+29,bx+24:bx+29) = diag( lambda(bc+20:bc+25).*  ( - ((1-x(bx+8:bx+13)).^2 + 2*ep^2)./ atemp ) );
        hess(bx+8:bx+13,bx+8:bx+13) = diag( lambda(bc+20:bc+25).* ( - (x(bx+24:bx+29).^2 + 2*ep^2) ./ atemp ) );
        % hess(bx+8:bx+13,bx+24:bx+29) = diag( lambda(bc+20:bc+25).* ( - x(bx+24:bx+29).*(1-x(bx+8:bx+13)) ./ atemp ) );
        hess(bx+24:bx+29,bx+8:bx+13) = diag( lambda(bc+20:bc+25).* ( - x(bx+24:bx+29).*(1-x(bx+8:bx+13)) ./ atemp ) );

        atemp = (x(bx+30:bx+35).^2 + x(bx+8:bx+13).^2 + 2*ep^2).^1.5;
        hess(bx+30:bx+35,bx+30:bx+35) = diag( lambda(bc+26:bc+31).* (- (x(bx+8:bx+13).^2 + 2*ep^2)./ atemp ) );
        hess(bx+8:bx+13,bx+8:bx+13) = hess(bx+8:bx+13,bx+8:bx+13) + diag( lambda(bc+26:bc+31).* (- (x(bx+30:bx+35).^2 + 2*ep^2)./ atemp));
        % hess(bx+8:bx+13,bx+30:bx+35) = diag( lambda(bc+26:bc+31).* ( x(bx+8:bx+13).*x(bx+30:bx+35) ./ atemp ) );
        hess(bx+30:bx+35,bx+8:bx+13) = diag( lambda(bc+26:bc+31).* ( x(bx+8:bx+13).*x(bx+30:bx+35) ./ atemp ) );

        bc = bc + 40;
        bx = bx + 39;
    end
end
end

function [hessstr] = hessian_structure_generation(dim_x,dim_nlcon,A_obj,A_control_P,ep,d,m,N)
x = 5 * rand(dim_x,1);
sigma = 5 * rand(1,1);
lambda = 5 * rand(dim_nlcon,1);
[hess] = hessian(x,sigma,lambda,dim_x,A_obj,A_control_P,ep,d,m,N);
hess(hess ~= 0) = 1;
hessstr = sparse(hess);
end

function [hessstr] = hessian_structure(hess)
hessstr = hess;
end
