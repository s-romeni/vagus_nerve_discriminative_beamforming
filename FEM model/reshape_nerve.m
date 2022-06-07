function circular_fascicles = reshape_nerve(circular_fascicles,h,R,delta)
%%-------------------------------------------------------------------------
% Function description: 
%%-------------------------------------------------------------------------
% Iteratevely move the fascicles to mimic the nerve re-organization after
% the implantation of a intraneural interface, in particular the TIME.
%%-------------------------------------------------------------------------
% Inputs: 
%%-------------------------------------------------------------------------
% • circular_fascicles [mm]: a numeric matrix containing the coordinates 
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
n_fasc = size(circular_fascicles,1);
fa = 1;
unf = 0;
while fa ~= 0
    intersection = zeros(n_fasc,4);
    interfasc = zeros(n_fasc,(n_fasc-1),4);
    for i = 1:n_fasc
        %%-----------------------------------------------------------------
        % Check if fascicle # i intercept the electrode
        [intersection(i,1:2),intersection(i,3:4)] = linecirc(0,h,circular_fascicles(i,1),circular_fascicles(i,2),circular_fascicles(i,3));
        %%-----------------------------------------------------------------
        % Move in the y axis each fascicle that intercept the electrode
        if not(isnan(sum(intersection(i,:))))
            s = sign(circular_fascicles(i,2)+h);
            if s == 0
                s = 1;
            end
            % do not move the fascicle outside the nerve
            if norm([abs(circular_fascicles(i,1)) abs(circular_fascicles(i,2)+s*delta*1e-2)]) < (R - circular_fascicles(i,3) - delta)
                circular_fascicles(i,2) = circular_fascicles(i,2)+s*delta*1e-2;
            end
        end
        %%-----------------------------------------------------------------
        % Check if the movement caused intercept other fascicle
        fasc = 1:n_fasc;
        fasc = setdiff(fasc,i);
        jj = 0;
        for j = 1:(n_fasc-1)
            jj = jj+1;
            % inter-fascicle interception
            [interfasc(i,jj,1:2),interfasc(i,jj,3:4)] = circcirc(circular_fascicles(i,1),circular_fascicles(i,2),circular_fascicles(i,3)+delta,circular_fascicles(fasc(j),1),circular_fascicles(fasc(j),2),circular_fascicles(fasc(j),3)+delta);
            
            if not(isnan(sum(interfasc(i,jj,:)))) % if there is an interception
                % control the position on the x axis
                s = sign(circular_fascicles(fasc(j),1)-circular_fascicles(i,1));
                % control the position on the y axis
                ss = sign(circular_fascicles(fasc(j),2)-circular_fascicles(i,2));

                % do not move the fascicle outside the nerve
                if norm([abs(circular_fascicles(fasc(j),1)+s*delta*1e-2) abs(circular_fascicles(fasc(j),2))]) < (R - circular_fascicles(fasc(j),3) - delta)
                    circular_fascicles(fasc(j),1) = circular_fascicles(fasc(j),1)+s*delta*1e-2;
                elseif norm([abs(circular_fascicles(fasc(j),1)) abs(circular_fascicles(fasc(j),2)+ss*delta*1e-2)]) < (R - circular_fascicles(fasc(j),3) - delta)
                    circular_fascicles(fasc(j),2) = circular_fascicles(fasc(j),2)+ss*delta*1e-2;
                else
                    unf = 1;
                    break
                end  
            end
        end
    end
    fa = sum(sum(not(isnan(intersection))));
    fa = fa+sum(sum(squeeze(sum(not(isnan(interfasc))))));
    if unf == 1
        fa = 0;
    end
end
