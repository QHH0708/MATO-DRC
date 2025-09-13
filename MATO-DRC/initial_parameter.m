%%% Initialization For Single Agent Trajectory Optimization with Obstacles

%%%  3 flat outputs, each integral chain is 3rd order
nu = [3;3;3];  %%% unchangeable
m = size(nu,1); % number of flat outputs


%%%   d=0 only position interior-point constraints(IPC), d=1 both position and velocity IPC
d = 0;


%%%   Weight matrix of the optimization objective
W = eye(3);


%%%   Number of polynomial trajectory segements
N = 6;


%%%  Time Instants
t0 = 0;  % start time
tN = 20; % end time

times = zeros(N+1,1);  % Split the time evenly
for n = 1:1:N
    times(n+1) = (tN/N) * n;
end


%%% Smoothing MPCC will be repeatedly solved multiple times with random initial points ...
%%% in order to calculate the average computational time 
try_times = 10;   


%%% the intial and end states
%%% [p1,v1,a1,p2,v2,a2,p3,v3,a3]^T
%%% p: position, v: velocity, a: acceleration
z0 = zeros(9,1);
z0(7) = 10;
zN = zeros(9,1);
zN(1) = 40; %%% x
zN(4) = 40; %%% y
zN(7) = 30; %%% z


%%% 2 ellipsoid obstacles and 2 zonotope obstacles
G_ellip = cell(2,1); %% generator matrices
c_ellip = cell(2,1); %% centers
G_zono = cell(2,1);
c_zono = cell(2,1);

G_ellip{1} = [6 0 0;0 6 0;0 0 18]; %% must be a 3*3 matrix
c_ellip{1} = [10;10;18];
G_ellip{2} = [6 0 0;0 6 0;0 0 18]; %% must be a 3*3 matrix
c_ellip{2} = [30;30;18];

G_zono{1} = [5 0 0;0 5 0;0 0 15]; %% must be a 3*3 matrix
c_zono{1} = [10;30;15];
G_zono{2} = [5 0 0;0 5 0;0 0 15]; %% must be a 3*3 matrix
c_zono{2} = [30;10;15];

% G_zono{1} = [5 0 7;0 5 0;0 -2 15]; %% must be a 3*3 matrix
% c_zono{1} = [10;30;15];
% G_zono{2} = [5 3 0;0 5 0;5 0 15]; %% must be a 3*3 matrix
% c_zono{2} = [30;10;15];



%%%  reachable set of the agent
%%%  in this simulation, only the size of agent is considered as the
%%%  reachable set, i.e., a ball withe the radius = 0.5
G_reach = cell(N,1);
c_reach = cell(N,1);
for n = 1:1:N
    G_reach{n} = 0.5*eye(3);
    c_reach{n} = zeros(3,1);
end
