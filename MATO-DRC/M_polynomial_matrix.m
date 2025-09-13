function [ret] = M_polynomial_matrix(nu,d,W,t)
%%%%%% nu vector,  weight matrix W^\prime = -1/2 inv(W)
W = - 0.5 * inv(W);

m = size(nu,1); % number of flat output
sum_nu = sum(nu);

WP = []; 
for i = 1:1:m
    raw_WP = [];
    for j = 1:1:m
        raw_WP = [raw_WP W(i,j)*polynomial_matrix(nu(j),nu(i),1,t)];
    end
    WP = [WP;raw_WP];
end

% disp(WP)

P_diag_1 = [];
for i = 1:1:m
    P_diag_1 = blkdiag(P_diag_1,polynomial_matrix(nu(i),0,1-nu(i),t));
end

P_diag_2 = [];
for i = 1:1:m
    P_diag_2 = blkdiag(P_diag_2,polynomial_matrix(nu(i),0,d-nu(i)+2,t));
end

ret = [P_diag_1 WP;zeros(sum_nu-m*(d+1),sum_nu) P_diag_2];

end