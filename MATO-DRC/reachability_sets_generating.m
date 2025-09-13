function [G_reach,c_reach,ReachableSet_c_lowdim,ReachableSet_G_lowdim] = reachability_sets_generating(tau,N,times)
%%%%% tau, tiem step; t_start t_end, start and end time

kv = 6;
kp = 8;


Am = [0 1;-kp -kv];
Bm = [0;1];

A = blkdiag(Am,Am,Am);
B = blkdiag(Bm,Bm,Bm);

c_E0 = zeros(6,1); %%%% The center of the initial state set
G_E0 = diag([1 0.5 1 0.5 1 0.5]); %%% The generator matrix of the initial state set; ellipsoid
r_E0 = max(eig(G_E0));
% plot_ellipsoid_3D(zeros(3,1),0.5*eye(3),'r',0.5,0.1,'-','r',0.5);

c_delta = zeros(3,1) + [0.1;0;0]; %%%% center of disturbance
G_delta = 0.5*eye(3); %%% generator matrix of disturbance; ellipsoid

c_Bdelta = B * c_delta;
G_Bdelta = B * G_delta;
r_Bdelta = sqrt(max(eig(G_Bdelta' * G_Bdelta)));

e_tauA = expm(tau*A);
a_tau = (exp(tau * norm(A,2)) - 1 - tau * norm(A,2)) * (r_E0 + r_Bdelta / norm(A,2));

c_ellip=cell(3,1);
G_ellip=cell(3,1);
c_ellip{1,1} = e_tauA * c_E0;
G_ellip{1,1} = e_tauA * G_E0;
c_ellip{2,1} = tau * c_Bdelta;
G_ellip{2,1} = tau * G_Bdelta;
c_ellip{3,1} = a_tau * zeros(6,1);
G_ellip{3,1} = a_tau * eye(6);
[c_temp,G_temp] = ellipsoid_MinkowskiSum_overapprox(3,c_ellip,G_ellip);

c_ellipnew=cell(2,1);
G_ellipnew=cell(2,1);
c_ellipnew{1,1} = c_E0;
G_ellipnew{1,1} = G_E0;
c_ellipnew{2,1} = c_temp;
G_ellipnew{2,1} = G_temp;
[c_Omega0,G_Omega0] = ellipsoid_ConvexHull(2,c_ellipnew,G_ellipnew,6);

Nt = (times(N+1) - times(1)) / tau;  %%%% number of discretized reachable sets
ReachableSet_c = cell(Nt,1);
ReachableSet_G = cell(Nt,1);
ReachableSet_c{1,1} = c_Omega0;
ReachableSet_G{1,1} = G_Omega0;

b_tau = (exp(tau * norm(A,2)) - 1 - tau * norm(A,2)) * (r_Bdelta / norm(A,2));
c_ellip{3,1} = b_tau * zeros(6,1);
G_ellip{3,1} = b_tau * eye(6);

c_Omegak = c_Omega0;
G_Omegak = G_Omega0;
for i = 2:1:Nt
    c_ellip{1,1} = e_tauA * c_Omegak;
    G_ellip{1,1} = e_tauA * G_Omegak;
    [c_Omegakk,G_Omegakk] = ellipsoid_MinkowskiSum_overapprox(3,c_ellip,G_ellip);
    ReachableSet_c{i,1} = c_Omegakk;
    ReachableSet_G{i,1} = G_Omegakk;
    c_Omegak = c_Omegakk;
    G_Omegak = G_Omegakk;
end

color = 'b';
alf = 0.2;
LineWidth = 0.1;
LineStyle = '-';
EdgeColor = 'b';
Edgealf = 0.2;

G_reach = cell(N,1);
c_reach = cell(N,1);

ReachableSet_c_lowdim = cell(Nt,1);
ReachableSet_G_lowdim = cell(Nt,1);

c_UAVsize = zeros(3,1);  %%%% size of a UAV is modelled as an ellipsoid
G_UAVsize = 0.5*eye(3);

M_selet_dim = blkdiag([1 0],[1 0],[1 0]);
bi = 0;
for n = 1:1:N
    % figure(n);
    Nn = (times(n+1)-times(n)) / tau;  %%% The discrete number of each trajectory segment
    for i = bi+1:1:bi+Nn
        c = M_selet_dim * ReachableSet_c{i,1};
        G = M_selet_dim * ReachableSet_G{i,1};  %%% Dimension reduction of ellipsoids
        [U,S,V] = svd(G);  %%%% Perform singular value decomposition on the generated matrix of ellipsoid after dimension reduction
        G_new = U * S(:,1:3);
        % %%%%% Add the size of a UAV to the reachable set
        [MinkSum_c,MinkSum_G] = ellipsoid_MinkowskiSum_overapprox(2,{c_UAVsize;c},{G_UAVsize;G_new});
        ReachableSet_c_lowdim{i,1} = MinkSum_c;
        ReachableSet_G_lowdim{i,1} = MinkSum_G;

        %%%%% Do not add the size of a UAV to the reachable set
        % ReachableSet_c_lowdim{i,1} = c;
        % ReachableSet_G_lowdim{i,1} = G_new;
    end
    [CovxH_c,CovxH_G] = ellipsoid_ConvexHull(Nn,ReachableSet_c_lowdim(bi+1:bi+Nn,1),ReachableSet_G_lowdim(bi+1:bi+Nn,1),3);
    c_reach{n,1} = CovxH_c;
    G_reach{n,1} = CovxH_G;
    % plot_ellipsoid_3D(CovxH_c,CovxH_G,'r',alf,LineWidth,LineStyle,'r',Edgealf);
    bi = bi + Nn;
end

% for i = 1:1:Nt
%     if mod(i,4) == 1
%     plot_ellipsoid_3D(ReachableSet_c_lowdim{i,1}+0.1*i,ReachableSet_G_lowdim{i,1},color,alf,LineWidth,LineStyle,EdgeColor,Edgealf); hold on
%     end
% end


% %%%%% test Minkowski Sum
% G1 = 2*eye(3);
% c1 = [0;0;0];
% G2 = [1 0 1;0 2 1;0 0 4];
% % G2 = 4*eye(1);
% c2 = [1;2;3];
% c_ellip=cell(2,1);
% G_ellip=cell(2,1);
% c_ellip{1,1} = c1;
% G_ellip{1,1} = G1;
% c_ellip{2,1} = c2;
% G_ellip{2,1} = G2;
% num_ellip = 2;

% [MinkSum_c,MinkSum_G] = ellipsoid_MinkowskiSum_overapprox(num_ellip,c_ellip,G_ellip);
% plot_ellipsoid_3D(c1,G1,'r',0.5,0.5,'-','r',0.5);
% plot_ellipsoid_3D(c2,G2,'b',0.5,0.5,'-','b',0.5);
% plot_ellipsoid_3D(MinkSum_c,MinkSum_G,'g',0.5,0.5,'-','g',0.5);
% %%%%% test ellipsoid convex hull
% [CovxH_c,CovxH_G] = ellipsoid_ConvexHull(num_ellip,c_ellip,G_ellip,3);
% plot_ellipsoid_3D(c1,G1,'r',0.5,0.5,'-','r',0.5);
% plot_ellipsoid_3D(c2,G2,'b',0.5,0.5,'-','b',0.5);
% plot_ellipsoid_3D(CovxH_c,CovxH_G,'g',0.5,0.5,'-','g',0.5);

end

function [MinkSum_c,MinkSum_G] = ellipsoid_MinkowskiSum_overapprox(num_ellip,c_ellip,G_ellip)
%%%% num_ellip, number of ellipsoids, c_ellip G_ellip, centers and generator matrices, cell(n,1)
MinkSum_c = zeros(size(c_ellip{1,1}));
for n = 1:1:num_ellip
    MinkSum_c = MinkSum_c + c_ellip{n,1};
end
temp1 = 0;
temp2 = zeros(size(G_ellip{1,1}));
for n = 1:1:num_ellip
    Q = G_ellip{n,1} * G_ellip{n,1}';
    sqrttrQ = sqrt(trace(Q));
    temp1 = temp1 + sqrttrQ;
    temp2 = temp2 + Q / sqrttrQ;
end
MinkSum_Q = temp1 * temp2;
MinkSum_G = MinkSum_Q^0.5;
end

function [CovxH_c,CovxH_G] = ellipsoid_ConvexHull(num_ellip,c_ellip,G_ellip,dim)
%%%% dim, dimension of an ellipsoid
yalmip('clear');
A0 = sdpvar(dim,dim,'symmetric');
b0 = sdpvar(dim,1);
t = sdpvar(num_ellip,1);
M0 = [A0 b0 zeros(dim,dim);b0' -1 b0';zeros(dim,dim) b0 -A0];
Constraints = [A0 >= 0.00001*eye(dim), t >= 0];
for n = 1:1:num_ellip
    An = inv(G_ellip{n,1} * G_ellip{n,1}');
    bn = - An * c_ellip{n,1};
    cn = c_ellip{n,1}' * An * c_ellip{n,1} - 1;
    Mn = [An bn zeros(dim,dim);bn' cn zeros(1,dim);zeros(dim,dim) zeros(dim,1) zeros(dim,dim)];
    Constraints = [Constraints; M0 - t(n) * Mn <= 0];
end
ObjFunction = - logdet(A0);
options = sdpsettings('verbose',1);  %%% verbose=0 no output information
sol = optimize(Constraints,ObjFunction,options);
A = value(A0);
b = value(b0);
CovxH_G = inv(A)^0.5;
CovxH_c = - A \ b;
end


