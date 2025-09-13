%%%%%%%%%%%%%% for multi agent with zonotope reachable set

t1 = datetime('now'); %%% for the computational time of zonotope reachable generation

nu = [3;3;3];
d = 0;
W = eye(3);
% times = [0;5;10;15;20];
N = 4;
m = size(nu,1);

t0 = 0;
tN = 20;

times = zeros(N+1,1);
for n = 1:1:N
    times(n+1) = (tN/N) * n;
end

Na = 5; %%% number of agents
tau = 0.05; %%% time step
% tau = 0.001;

try_times = 2;


z0 = cell(1,Na);
zN = cell(1,Na);
for a = 1:1:Na
    z0{1,a} = zeros(9,1);
    zN{1,a} = zeros(9,1);
end


points = [0 35 15;
          40 20 30;
          0 5 10;
          23 40 30;
          27 0 20;
          0 35 15];

for a = 1:1:Na
    for i = 1:1:3
        z0{1,a}((i-1)*3+1) = points(a,i);
        zN{1,a}((i-1)*3+1) = points(a+1,i);
    end
end

G_ellip = cell(2,1);
c_ellip = cell(2,1);
G_zono = cell(2,1);
c_zono = cell(2,1);

G_ellip{1} = [6 0 0;0 6 0;0 0 18];
c_ellip{1} = [10;10;18];
G_ellip{2} = [6 0 0;0 6 0;0 0 18];
c_ellip{2} = [30;30;18];

G_zono{1} = [5 0 0;0 5 0;0 0 15];
c_zono{1} = [10;30;15];
G_zono{2} = [5 0 0;0 5 0;0 0 15];
c_zono{2} = [30;10;15];

G_reach = cell(N,Na);
c_reach = cell(N,Na);

[G_temp,c_temp,Reachset_c,Reachset_G] = reachability_sets_generating_zonotope(tau,N,times);

% for n = 1:1:N
%     disp(G_temp{n,1});
%     disp(c_temp{n,1});
% end
for a = 1:1:Na  %%%% each agent's reachable set is the same in identical time interval
    for n = 1:1:N
        G_reach{n,a} = G_temp{n,1};
        c_reach{n,a} = c_temp{n,1};
    end
end

t2 = datetime('now');
computing_time_zonoRS_generation = seconds(t2 - t1);
disp(computing_time_zonoRS_generation);