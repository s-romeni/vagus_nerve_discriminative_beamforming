function data = physiomodel(t,xyz,A,fresp,rr,dt,thetai,a,b)
%%-------------------------------------------------------------------------
% Function description: 
%%-------------------------------------------------------------------------
% This function generate synthetic BP (or eventually ECG and RESP) given the model parameters.
%%-------------------------------------------------------------------------
% Inputs: 
%%-------------------------------------------------------------------------
% • t: time.
% • xyz: initial conditions.
% • A: desired signal amplitude.
% • fresp: respiratory frequency.
% • rr: inter-beat interval process.
% • dt: time resolution.
% • thetai: waves phase parameters.
% • a: waves amplitude parameters.
% • b: waves width parameters.
%%-------------------------------------------------------------------------
% Outputs:
%%-------------------------------------------------------------------------
% • data: synthetic physiological data.
%%-------------------------------------------------------------------------
% Source: 
%%-------------------------------------------------------------------------
% Source: https://archive.physionet.org/physiotools/ecgsyn/Matlab/
%%-------------------------------------------------------------------------
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

