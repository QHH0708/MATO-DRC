function [s] = plot_ellipsoid_3D(c,G,color,alf,LineWidth,LineStyle,EdgeColor,Edgealf,accuracy)

r = 1;

vert_cir=[];
% n2=12;  %Divide the angle into equal parts
% n1=6;
n2 = accuracy;
n1 = n2 /2;
for i=0:1:n1
    phi = i * pi/n1;
    rz = r * cos(phi);
    r2 = r * sin(phi);
    for j = 1:1:n2
        theta = j * 2*pi/n2;
        vert_cir=[vert_cir;r2*cos(theta) r2*sin(theta) rz];
    end
end
set_cir=G*Polyhedron(vert_cir)+c;

s=plot(set_cir,'color',color);  hold on

alpha(s,alf);
set(s,'LineWidth',LineWidth,'LineStyle',LineStyle,'EdgeColor',EdgeColor,'EdgeAlpha',Edgealf);


end