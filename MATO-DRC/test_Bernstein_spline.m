function test_Bernstein_spline()

n = 6;

l = 0;
u = 4;

M1 = zeros(n+1,n+1);

for i = 0:1:n
    c = nchoosek(n,i);

    for j = n-i:-1:0
        c1 = nchoosek(n-i,j);
        m = c * c1 * (-1)^j;
        M1(i+1,n+1-j-i) = m;
    end

end

disp(M1);

% A = [-3 7 -20 4;
%      1 -5 10 25;
%      2 2 -8 -5];

% A = [0 7 -20 4;
%      0 -5 10 25;
%      0 2 -8 -5];

% A = [7 -20 4;
%      -5 10 25;
%      2 -8 -5];

A = [0 0 7 -20 4 2 -3;
     0 0 -5 10 25 -1 5;
     0 0 2 -8 -5 1 2];

A_new = zeros(size(A));

for k = 1:1:3
    for i = n:-1:0
        for j = i:-1:0
            c = nchoosek(i,j);
            A_new(k,n-j+1) = A_new(k,n-j+1) + A(k,n-i+1) * c * l^(i-j) * (u-l)^(j);
        end
    end
end

disp(A_new);

control_P = M1 \ A_new';
disp(control_P)

np = (u-l)*10 + 1;
points = zeros(3,np);

for i=1:1:np
    points(1,i) = polyval(A(1,:),0.1*(i-1));
    points(2,i) = polyval(A(2,:),0.1*(i-1));
    points(3,i) = polyval(A(3,:),0.1*(i-1));
end

plot3(points(1,:),points(2,:),points(3,:)); hold on

control_P = control_P';
disp(control_P)
plot3(control_P(1,:),control_P(2,:),control_P(3,:),'-o')

end

