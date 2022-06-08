% This function has been used to generate synthetic BP
% It can be used also for ECG and RESP
%%-------------------------------------------------------------------------
function data = physiomodel(t,xyz,A,fresp,rr,dt,thetai,a,b)
theta = atan2(xyz(2),xyz(1));
alpha = 1 - sqrt(xyz(1)^2+xyz(2)^2);
i = t/dt;
omega = 2*pi/rr(floor(i));
z0 = A*sin(2*pi*fresp*t);
dx = alpha*xyz(1) - omega*xyz(2);
dy = alpha*xyz(2) + omega*xyz(1);
deltatheta = rem((theta - thetai),2*pi);
dz = - sum(a.*deltatheta.*exp(-0.5*((deltatheta./b).^2))) - (xyz(3)-z0);
data = [dx;dy;dz];
end

