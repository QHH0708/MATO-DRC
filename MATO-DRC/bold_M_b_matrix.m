function [ret_M,ret_b] = bold_M_b_matrix(nu,d,W,times,N,z0,zN,zd)
%%%%%% nu vector;  times: time instant vector
% N = size(times,1) - 1;
%%% z0 zN start and end of a trajectory
%%% zd: (N-1) row  m*(d+1) col matrix, path point constraints

m = size(nu,1);  %%% number of flat outputs

sum_nu = sum(nu);

ret_M = zeros(2*N*sum_nu,2*N*sum_nu);

R0 = R_polynomial_matrix(nu,W,times(1));
RN = R_polynomial_matrix(nu,W,times(N+1));

ret_M(1:sum_nu,1:2*sum_nu) = R0;
ret_M((2*N-1)*sum_nu+1:2*N*sum_nu,(2*N-2)*sum_nu+1:2*N*sum_nu) = RN;

raw_base = sum(nu);
col_base = 0;
for i = 1:1:N-1
    Qi = Q_polynomial_matrix(nu,d,W,times(i+1));
    Mi = M_polynomial_matrix(nu,d,W,times(i+1));
    ret_M(raw_base+1:raw_base+m*(d+1),col_base+1:col_base+2*sum_nu) = Qi;
    ret_M(raw_base+m*(d+1)+1:raw_base+2*sum_nu,col_base+1:col_base+4*sum_nu) = [Mi -Mi]; 
    raw_base = raw_base + 2*sum_nu;
    col_base = col_base + 2*sum_nu;
end

ret_b = zeros(2*N*sum_nu,1);

ret_b(1:sum_nu,1) = z0;
ret_b((2*N-1)*sum_nu+1:2*N*sum_nu,1) = zN;

raw_base1 = sum_nu;
for i = 1:1:N-1
    ret_b(raw_base1+1:raw_base1+m*(d+1),1) = zd(i,:)';
    raw_base1 = raw_base1 + 2*sum_nu;
end

end