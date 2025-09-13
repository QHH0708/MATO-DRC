function test_Ts(bar_nu_z,t_start,t_end)
Ts = Ts_time_scaling_matrix_Bernstein(bar_nu_z,t_start,t_end);
test_Ts_result = [];
n = (t_end - t_start) / 0.1 + 1;
for i = 1:1:n
    t = t_start + (i-1)*0.1;
    t_prime = (t - t_start) / (t_end - t_start);
    test_Ts_result = [test_Ts_result;polynomial_vector(bar_nu_z-1,t)-polynomial_vector(bar_nu_z-1,t_prime)*Ts];
end
disp(test_Ts_result);
end

function [ret] = polynomial_vector(order,t)
%%%%%% order = hightest power - 1
ret = zeros(1, order + 1);
for i = 0:1:order
    ret(1,i+1) = t^i / factorial(i);
end
end