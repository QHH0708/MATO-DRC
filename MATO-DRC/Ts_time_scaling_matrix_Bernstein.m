function [Ts] = Ts_time_scaling_matrix_Bernstein(bar_nu_z,t_start,t_end)
%%% polynomial vector P(t) = P(t^\prime) Ts

Ts = zeros(bar_nu_z,bar_nu_z);

for i = 0:1:bar_nu_z - 1  %%% i for col of Ts
    fi = factorial(i);
    for j = 0:1:i
        Ts(j+1,i+1) = nchoosek(i,j) * (t_end - t_start)^j * t_start^(i-j) * factorial(j) / fi;
    end
end

end

