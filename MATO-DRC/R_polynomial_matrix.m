function [ret] = R_polynomial_matrix(nu,W,t)
%%%%%% nu vector， W: weight matrix, W^\prime = -1/2 W
W = - 0.5 * inv(W);

m = size(nu,1); % number of flat output
% sum_nu = sum(nu);

WP = []; 
for i = 1:1:m
    raw_WP = [];
    for j = 1:1:m
        raw_WP = [raw_WP W(i,j)*polynomial_matrix(nu(j),nu(i),1,t)];
    end
    WP = [WP;raw_WP];
end

% disp(WP)

P_diag = [];
for i = 1:1:m
    P_diag = blkdiag(P_diag,polynomial_matrix(nu(i),0,1-nu(i),t));
end

ret = [P_diag WP];

end