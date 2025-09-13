function [ret] = polynomial_matrix(original_order,s1,s2,t) 
%%%% P_(original_order)^[s1,s2] (t)
%%%% original_order = highest power + 1
%%%%%%% s1 >= s2
%    %%%%%% order = highest power - 1

% [ret] = polynomial_vector(order,t);

ret = zeros(s1-s2+1,original_order);
for i = 1:1:s1-s2+1
    if s1 - i + 1 >= 0
        ret(i,:) = polynomial_vector_integral(original_order,s1-i+1,t);
    else
        ret(i,:) = polynomial_vector_derivative(original_order,-(s1-i+1),t);
    end
end
end

% function [ret] = polynomial_vector(order,t)
% %%%%%% order = highest power - 1
% ret = zeros(1, order + 1);
% for i = 0:1:order
%     ret(1,i+1) = t^i / factorial(i);
% end
% end

function [ret] = polynomial_vector_integral(original_order,integral_order,t)
%%%%%%% origianl_order, original order of polynomial; integral_order, additional integral order
ret = zeros(1,original_order);
for i = 1:1:original_order
    ret(1,i) = t^(i - 1 + integral_order) / factorial(i - 1 + integral_order);
end
end

function [ret] = polynomial_vector_derivative(original_order,deriv_order,t)
%%%%%%%% origianl_order, original order of polynomial;  deriv_order,  additional derivative order
ret = zeros(1,original_order);
for i = deriv_order+1:1:original_order
    ret(1,i) = t^(i - 1 - deriv_order) / factorial(i - 1 - deriv_order);
end
end


