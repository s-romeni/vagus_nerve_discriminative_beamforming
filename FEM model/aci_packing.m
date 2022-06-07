function [centers,r] = aci_packing(R, rmax, rmin, delta, epsilon, N)
%%-------------------------------------------------------------------------
% Function description:  
%%-------------------------------------------------------------------------
% A-priori Check for Intersections packing of circular objects
% in a circular region.
%%-------------------------------------------------------------------------
% Inputs:
%%-------------------------------------------------------------------------
% • R: radius of the region to be packed;
% • rmax: maximum objects radius
% • rmin: minimum objects radius
% • delta: minimum spacing between packed objects;
% • epsilon: number of grid refinements.
%%-------------------------------------------------------------------------
% Outputs
% • centers: Nx2 array with object center positions.
% • r: Nx1 array of radii of objects.
%%-------------------------------------------------------------------------
% Author: 
%%-------------------------------------------------------------------------
% Simone Romeni @TNE, EPFL
%%-------------------------------------------------------------------------
r = (rmax-rmin)*rand(N,1) + rmin;
xmin = -R; xmax = R; ymin = -R; ymax = R;
x = xmin:epsilon:xmax;
y = ymin:epsilon:ymax;
[x,y] = meshgrid(x,y);
% count = 0;
% maxcount = 10;
iter = 1;
niter = 1;
while iter == 1
    centers = zeros(N,2);
    for i = 1:N
        sel = x.^2 + y.^2 < (R-r(i)-delta)^2;
        for j = 1:i-1
            sel = sel & ((x-centers(j,1)).^2 + (y-centers(j,2)).^2 > (r(j)+r(i)+delta)^2);
        end
        
        % Try packing maxcount times
        %     if ~sum(sum(selaux))
        %         warning('Packing interrupted!')
        %         count = count + 1;
        %         if count >= maxcount
        %             error('No try went well...')
        %         end
        %     end
        
        [row,col] = find(sel);
        if isempty(row) || isempty(col)
            niter = niter+1;
            if niter > 300
                r = (rmax-rmin)*rand(N,1) + rmin;
                niter = 1;
            end
            break
        else
            nrand = randi(length(row));
            centers(i,:) = [x(row(nrand),col(nrand)),y(row(nrand),col(nrand))];
            if i == N
                iter = 0;
            end
        end
    end
end
end
