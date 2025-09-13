function [ret] = R_reorganized_trj_parameter(nu,N,bold_M,bold_b,W)
%%%%%% nu vector

sum_nu = sum(nu);
m = size(nu,1);

bold_c_all = bold_M \ bold_b;

[reorganize_matrices] = reorganize_each_flatoutput(nu,W);

for n = 1:1:N   %%% each piece of trajectory
    bold_c_n = bold_c_all((n-1)*2*sum_nu+1:n*2*sum_nu);
    
end

end

function [ret] = reorganize_each_flatoutput(nu,W)   
%%%% compute reorganize matrix for each flat output
%%%% ret{i,1} * old_parameter = new_parameter
%%%%%% nu vector

m = size(nu,1);
sum_nu = sum(nu);

bar_nu_z = 2 * max(nu);

ret = cell(m,1);  %%% reorganize matrix of each flat output

raw_base_i = 0;
for i = 1:1:m
    reorganize_matrix_i = zeros(bar_nu_z,2*sum_nu);

    reorganize_matrix_i(1:nu(i),raw_base_i+1:raw_base_i+nu(i)) = eye(nu(i));
    raw_base_i = raw_base_i + nu(i);

    col_base_j = sum_nu;
    for j = 1:1:m
        temp = reorganize_matrix_i(nu(i)+1:nu(i)+nu(j),col_base_j+1:col_base_j+nu(j));
        reorganize_matrix_i(nu(i)+1:nu(i)+nu(j),col_base_j+1:col_base_j+nu(j)) = temp + W(i,j)*eye(nu(j));
        col_base_j = col_base_j + nu(j);
    end

    ret{i,1} = reorganize_matrix_i;
end

end