function [Wint] = Wint_matrix_ofobjectivefunction(nu,W,times,N)
%%% W weight matrix

W = inv(W);  %%% inverse

m = size(nu,1);
sum_nu = sum(nu);

Wint = zeros(2*N*sum_nu,2*N*sum_nu);

base_n = sum_nu;
for n = 1:1:N
    t_start = times(n);
    t_end = times(n+1);
    Wint_n = zeros(sum_nu,sum_nu);
    base_i = 0;
    for i = 1:1:m
        base_j = 0;
        for j = 1:1:m
            Wint_n(base_i+1:base_i+nu(i),base_j+1:base_j+nu(j))=W(i,j)*(integral_PTP(nu(i),nu(j),t_end)-integral_PTP(nu(i),nu(j),t_start));
            base_j = base_j + nu(j);
        end
        base_i = base_i + nu(i);
    end
    Wint(base_n+1:base_n+sum_nu,base_n+1:base_n+sum_nu) = Wint_n;
    base_n = base_n + 2*sum_nu;
end

Wint = Wint / 4; 

end



function [ret] = integral_PTP(nu_i,nu_j,t)
%%%% Multiply polynomial vectors and integrate once
ret = zeros(nu_i,nu_j);
for i = 1:1:nu_i
    for j = 1:1:nu_j
        ret(i,j) = t^(i+j-1) / (factorial(i-1) * factorial(j-1) * (i+j-1));
    end
end
end