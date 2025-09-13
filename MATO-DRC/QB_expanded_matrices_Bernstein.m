function [QB] = QB_expanded_matrices_Bernstein(bar_nu_z)
%%% bar_nu_z maximal order+1 of the trajectory
%%% Q_B Polynomial_vector ^ T
%%% the order of Q_B col = 0 1 2 ... bar_nu_z-1  

QB = zeros(bar_nu_z,bar_nu_z);

for i = 0:1:bar_nu_z - 1
    c = nchoosek(bar_nu_z - 1,i);
    for j = bar_nu_z - 1 - i:-1:0
        c1 = nchoosek(bar_nu_z - 1 - i,j);
        m = c * c1 * (-1)^j;
        QB(i+1,i+j+1) = m * factorial(i+j);
    end
end


end