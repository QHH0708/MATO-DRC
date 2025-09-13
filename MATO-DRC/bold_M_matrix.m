function [ret,zd_matrix_bold_b,vector_z0N] = bold_M_matrix(nu,d,W,times,N,z0,zN)
%%%%%% nu vector,  times, time instant vector
% N = size(times,1) - 1;

m = size(nu,1);  %%% number of flat outputs

sum_nu = sum(nu);

ret = zeros(2*N*sum_nu,2*N*sum_nu);

R0 = R_polynomial_matrix(nu,W,times(1));
RN = R_polynomial_matrix(nu,W,times(N+1));

ret(1:sum_nu,1:2*sum_nu) = R0;
ret((2*N-1)*sum_nu+1:2*N*sum_nu,(2*N-2)*sum_nu+1:2*N*sum_nu) = RN;

raw_base = sum(nu);
col_base = 0;
for i = 1:1:N-1
    Qi = Q_polynomial_matrix(nu,d,W,times(i+1));
    Mi = M_polynomial_matrix(nu,d,W,times(i+1));
    ret(raw_base+1:raw_base+m*(d+1),col_base+1:col_base+2*sum_nu) = Qi;
    ret(raw_base+m*(d+1)+1:raw_base+2*sum_nu,col_base+1:col_base+4*sum_nu) = [Mi -Mi]; 
    raw_base = raw_base + 2*sum_nu;
    col_base = col_base + 2*sum_nu;
end

zd_matrix_bold_b = zeros(2*N*sum_nu,(N-1)*m*(d+1));
raw_base = sum_nu;
for i = 1:1:N-1
    zd_matrix_bold_b(raw_base+1:raw_base+m*(d+1),(i-1)*m*(d+1)+1:i*m*(d+1)) = eye(m*(d+1));
    raw_base = raw_base + 2*sum_nu;
end

vector_z0N = zeros(2*N*sum_nu,1);
vector_z0N(1:sum_nu,1) = z0;
vector_z0N((2*N-1)*sum_nu+1:2*N*sum_nu,1) =  zN;

end