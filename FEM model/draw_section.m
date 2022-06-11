function draw_section(R,circular_fascicles_TIME,circular_fascicles,l_shaft,pos,d_as,h_as)
%%-------------------------------------------------------------------------
% Function description:
%%-------------------------------------------------------------------------
% Draw the nerve section to check the topography and the reorganization
%%-------------------------------------------------------------------------
% Inputs:
%%-------------------------------------------------------------------------
% • circular_fascicles [mm]: a numeric matrix containing the coordinate
% of each fascicle in a bidimensional plane [x,y] and his radius.
% • h [mm]: minimum height for fascicles boundary.
% • R [mm]: nerve radius.
% • delta [mm]: minimum spacing between electrode and fascicles boundary.
%%-------------------------------------------------------------------------
% Outputs:
%%-------------------------------------------------------------------------
% • circular_fascicles [mm]: a numeric matrix containing the (updated)
% coordinates of each fascicle in a bidimensional plane [x,y] and his
% radius.
%%-------------------------------------------------------------------------
% Authors: 
%%-------------------------------------------------------------------------
% Andrea Pitzus @TNE, SSSA // @MeDSP, UniCa & Simone Romeni @TNE, EPFL
%%-------------------------------------------------------------------------
fasc_number = size(circular_fascicles, 1);
n = fasc_number/2;
theta = 0:0.0105:2*pi+0.0105;
figure
hold on
for i = 1:fasc_number
    xi = (circular_fascicles(i, 3)*cos(theta)+circular_fascicles(i, 1));
    yi = (circular_fascicles(i, 3)*sin(theta)+circular_fascicles(i, 2));
    plot(xi,yi,'-k','LineWidth',2)
end
for i = 1:fasc_number
    xi = (circular_fascicles_TIME(i, 3)*cos(theta)+circular_fascicles_TIME(i, 1));
    yi = (circular_fascicles_TIME(i, 3)*sin(theta)+circular_fascicles_TIME(i, 2));
    plot(xi,yi,'-b','LineWidth',2)
end
xi = R*cos(theta);
yi = R*sin(theta);
plot(xi,yi,'-k','LineWidth',3);
for i = 1:2*n
    if i <= n
        xi = pos(i,1);
        yi = pos(i,2);
        rectangle('Position',[xi-d_as yi-2*h_as d_as h_as])
    elseif i > n
        xi = pos(i,1);
        yi = pos(i,2);
        rectangle('Position',[xi-d_as yi-h_as d_as h_as])
    end
end
rectangle('Position',[-R+l_cc/2 (-h_shaft/2-h_as) l_shaft h_shaft])
xlabel('x axis (mm)')
ylabel('y axis (mm)')
xticks(1e-3*[-1 -0.5 0 0.5 1])
xticklabels({'-1','-0.5','0','0.5','1'})
yticks(1e-3*[-1 -0.5 0 0.5 1])
yticklabels({'-1','-0.5','0','0.5','1'})
title('Human-like Vagus Section')
set(gca,'FontName','SansSerif')
set(gca,'FontSize',17)
axis(gca,'equal')
box on
end
