function [G_reach,c_reach,ReachableSet_c_lowdim,ReachableSet_G_lowdim] = reachability_sets_generating_zonotope(tau,N,times)
%%%%% tau, time step;t_start t_end, start and end time

kv = 6;
kp = 8;

Am = [0 1;-kp -kv];
Bm = [0;1];

A = blkdiag(Am,Am,Am);
B = blkdiag(Bm,Bm,Bm);

c_E0 = zeros(6,1); %%%% The center of the initial state set
G_E0 = diag([1 0.5 1 0.5 1 0.5]); %%% The generator matrix of the initial state set; zonotope
r_E0 = 1;

% plot_constrained_zonotope_3D(c,G,Ae,be,color,alf,LineWidth,LineStyle,EdgeColor,Edgealf)

c_delta = zeros(3,1) + [0.1;0;0]; %%%% center of disturbance
G_delta = 0.5*eye(3); %%% generator matrix of disturbance; zonotope
% G_delta = [0.4 0 0.1;
%           0.1 0.4 0;
%           0.1 0.3 0.1;];

c_Bdelta = B * c_delta;
G_Bdelta = B * G_delta;
r_Bdelta = 0.5;

e_tauA = expm(tau*A);
a_tau = (exp(tau * norm(A,inf)) - 1 - tau * norm(A,inf)) * (r_E0 + r_Bdelta / norm(A,inf));

c_z=cell(3,1);
G_z=cell(3,1);
c_z{1,1} = e_tauA * c_E0;
G_z{1,1} = e_tauA * G_E0;
c_z{2,1} = tau * c_Bdelta;
G_z{2,1} = tau * G_Bdelta;
c_z{3,1} = a_tau * zeros(6,1);
G_z{3,1} = a_tau * eye(6);
[c_temp,G_temp] = zonotope_MinkowskiSum(3,c_z,G_z);

c_znew=cell(2,1);
G_znew=cell(2,1);
c_znew{1,1} = c_E0;
G_znew{1,1} = G_E0;
c_znew{2,1} = c_temp;
G_znew{2,1} = G_temp;
[c_Omega0,G_Omega0] = zonotope_ConvexHull(2,c_znew,G_znew,6);

Nt = (times(N+1) - times(1)) / tau;  %%%% number of discretized reachable sets
ReachableSet_c = cell(Nt,1);
ReachableSet_G = cell(Nt,1);
ReachableSet_c{1,1} = c_Omega0;
ReachableSet_G{1,1} = G_Omega0;

b_tau = (exp(tau * norm(A,inf)) - 1 - tau * norm(A,inf)) * (r_Bdelta / norm(A,inf));
c_z{3,1} = b_tau * zeros(6,1);
G_z{3,1} = b_tau * eye(6);

c_Omegak = c_Omega0;
G_Omegak = G_Omega0;
for i = 2:1:Nt
    c_z{1,1} = e_tauA * c_Omegak;
    G_z{1,1} = e_tauA * G_Omegak;
    [c_Omegakk,G_Omegakk] = zonotope_MinkowskiSum(3,c_z,G_z);

    if mod(i,10) == 1  %%% order reduction of zonotope for every 10 steps
        G_Omegakk = zono_order_reduction(G_Omegakk, 3); %%% reduced to 3rd order
    end

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

ReachableSet_c_lowdim = cell(Nt,1);  %%% reachable sets in position space
ReachableSet_G_lowdim = cell(Nt,1);

% a = pi/4;
% RM_x = [1 0 0;
%         0 cos(a) -sin(a);
%         0 sin(a) cos(a)];
% RM_y = [cos(a) 0 sin(a);
%         0 1 0;
%         -sin(a) 0 cos(a)];
% a = pi/2;
% RM_z = [cos(a) -sin(a) 0;
%         sin(a) cos(a) 0;
%         0 0 1];
c_UAVsize = zeros(3,1);  %%%% size of a UAV is modelled as a zonotope
G_UAVsize = 0.5*eye(3);
% disp(G_UAVsize);
% plot_constrained_zonotope_3D(c_UAVsize,G_UAVsize,[],[],'r',alf,LineWidth,LineStyle,EdgeColor,Edgealf);

M_selet_dim = blkdiag([1 0],[1 0],[1 0]);
bi = 0;
for n = 1:1:N
    % figure(n);
    Nn = (times(n+1)-times(n)) / tau;  %%% The discrete number of each trajectory segment
    for i = bi+1:1:bi+Nn
        c = M_selet_dim * ReachableSet_c{i,1};
        G = M_selet_dim * ReachableSet_G{i,1};  

        % %%%%% add the size of quadrotor to reachable sets
        [MinkSum_c,MinkSum_G] = zonotope_MinkowskiSum(2,{c_UAVsize;c},{G_UAVsize;G});
        ReachableSet_c_lowdim{i,1} = MinkSum_c;
        ReachableSet_G_lowdim{i,1} = MinkSum_G;

    end
    [CovxH_c,CovxH_G] = zonotope_ConvexHull(Nn,ReachableSet_c_lowdim(bi+1:bi+Nn,1),ReachableSet_G_lowdim(bi+1:bi+Nn,1),3);
    c_reach{n,1} = CovxH_c;
    G_reach{n,1} = CovxH_G;

    % plot_ellipsoid_3D(CovxH_c,CovxH_G,'r',alf,LineWidth,LineStyle,'r',Edgealf);
    % plot_constrained_zonotope_3D(CovxH_c,CovxH_G,[],[],color,alf,LineWidth,LineStyle,EdgeColor,Edgealf)

    bi = bi + Nn;
end

%%% order reduction of the position reachable sets, accelerate the generation of figures
for i=1:1:Nt
    ReachableSet_G_lowdim{i,1} = zono_order_reduction(ReachableSet_G_lowdim{i,1}, 2); 
end


end


function [MinkSum_c, MinkSum_G] = zonotope_MinkowskiSum(num_zono,c_zono,G_zono)

MinkSum_c = zeros(size(c_zono{1}));
MinkSum_G = [];
for i = 1:1:num_zono
    MinkSum_c = MinkSum_c + c_zono{i};
    MinkSum_G = [MinkSum_G G_zono{i}];
end

end

function [CovxH_c, CovxH_G] = zonotope_ConvexHull(num_zono,c_zono,G_zono,dim)
%%% this convex hull is over-approximated

low_bound = zeros(num_zono, dim);
upper_bound = zeros(num_zono, dim);

for i = 1:1:num_zono
    for d = 1:1:dim
        interval_d = norm(G_zono{i}(d, :), 1);
        low_bound(i, d) = c_zono{i}(d, 1) - interval_d;
        upper_bound(i, d) = c_zono{i}(d, 1) + interval_d;
    end
end

low_bound_ConvexHull = min(low_bound, [], 1);
upper_bound_ConvexHull = max(upper_bound, [], 1);

CovxH_c = (upper_bound_ConvexHull + low_bound_ConvexHull)' / 2;
CovxH_G = diag((upper_bound_ConvexHull - low_bound_ConvexHull) / 2);

end


function [Gnew] = zono_order_reduction(G,order) %input: G, generator matrix of zonotope, order, the order after reduced
    row=size(G,1); 
    if order==1
        R=zeros(row,1);
        for i=1:1:row
            R(i)=norm(G(i,:),1);
        end
        Gnew=diag(R);
    end
    
    if order>1
        colnorm=zeros(1,size(G,2));
        for i=1:1:size(G,2)
            colnorm(i)=norm(G(:,i),2);
        end
        for i=1:1:size(G,2)-1
            for j=i+1:1:size(G,2)
                if colnorm(j)>colnorm(i)
                    colt=colnorm(j);
                    colnorm(j)=colnorm(i);
                    colnorm(i)=colt;
                    
                    ctemp=G(:,j);
                    G(:,j)=G(:,i);
                    G(:,i)=ctemp;
                end
            end
        end
        
        G1=G(:,1:order*row-row);
        R1=zeros(row,1);
        R2=zeros(row,1);
        for i=1:1:row
            R1(i)=norm(G1(i,:),1);
            R2(i)=norm(G(i,:),1);
        end
        Gnew=[G1 diag(R2-R1)];
    end
end
