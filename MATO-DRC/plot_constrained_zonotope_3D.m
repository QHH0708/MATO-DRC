function [set_con_zono]=plot_constrained_zonotope_3D(c,G,Ae,be,color,alf,LineWidth,LineStyle,EdgeColor,Edgealf)
%%% alf: transparency
[~,n]=size(G);
A=[eye(n);-eye(n)];
b=[ones(n,1);ones(n,1)];
box=Polyhedron('A',A,'b',b,'Ae',Ae,'be',be);
set_con_zono=c+G*box;

s=plot(set_con_zono,'color',color);  hold on
alpha(s,alf);

set(s,'LineWidth',LineWidth,'LineStyle',LineStyle,'EdgeColor',EdgeColor,'EdgeAlpha',Edgealf);




end