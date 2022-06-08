clc
clear
close all
%%-------------------------------------------------------------------------
% General info: generate synthetic blood pressure signal using a realistic 
% coupled nonlinear artificial ECG, BP and respiratory signal generator
% for assessing noise performance of biomedical signal processing
% algorithms
%%-------------------------------------------------------------------------
% The scripts is the implementation of the following bibliographic source:
% G. D. Clifford and P. E. McSharry, "A realistic coupled nonlinear 
% artificial ECG, BP, and respiratory signal generator for assessing noise 
% performance of biomedical signal processing algorithms," Maspalomas, 
% Gran Canaria Island, Spain, May 2004, p. 290. doi: 10.1117/12.544525.
%%-------------------------------------------------------------------------
% Source: https://archive.physionet.org/physiotools/ecgsyn/Matlab/
%%-------------------------------------------------------------------------
% Authors: 
%%-------------------------------------------------------------------------
% Andrea Pitzus @TNE, SSSA // @MeDSP, UniCa & Simone Romeni @TNE, EPFL
%%-------------------------------------------------------------------------

%%-------------------------------------------------------------------------
% BP
%%-------------------------------------------------------------------------
% index = [PBP QBP RBP SBP TBP]
thetai = [-(5/12)*pi -(1/36)*pi 0 (1/18)*pi (4/9)*pi];
a = [0 0 0.45 0.25 0.45];
b = [0.25 0.1 0.3 0.5 0.3];
%%-------------------------------------------------------------------------
% Runge Kutta 4th
%%-------------------------------------------------------------------------
dt = 0.025*1e-3;
t = (dt:dt:41.25)';
A = 0.15*1e-3; % V
fresp = 0.2667; % 16 resp cycle / min
seed = 1;
rng(seed,'twister');
%%-------------------------------------------------------------------------
% Define parameters for RR process 
%%-------------------------------------------------------------------------
% flo and fhi correspond to the Mayer waves
% and respiratory rate respectively
%%-------------------------------------------------------------------------
flo = 0.1;
fhi = 0.2667;
flostd = 0.001;
fhistd = 0.01;
lfhfratio = 0.5;
%%-------------------------------------------------------------------------
hrmean = 75; % bpm
hrstd = 2;
rrmean = (60/hrmean);
sampfreqrr = 1;
trr = 1/sampfreqrr; 
N = floor(2*hrmean/3); % approximate number of heartbeats
Nrr = 2^(ceil(log2(N*rrmean/trr)));
%%-------------------------------------------------------------------------
rr = rrprocess(flo,fhi,flostd,fhistd,lfhfratio,hrmean,hrstd,sampfreqrr,Nrr);
rr = interp(rr,1/dt);
%%-------------------------------------------------------------------------
xyz = [1,0,0.013];
[t,data] = ode45('physiomodel',t,xyz,[],A,fresp,rr,dt,thetai,a,b);
bp = data(:,3);
ot = (dt:dt:40)';
difft = length(t)-length(ot);
t = t(1:length(ot));
bp = bp(difft+1:end);
save('syn_bp.mat','bp','t')
%%-------------------------------------------------------------------------
% ECG
%%-------------------------------------------------------------------------
thetai = [-1.22173 -0.261799 0 0.261799 1.74533];
a = [1.2 -5 30 -7.5 0.75];
b = [0.25 0.1 0.1 0.1 0.4];
%%-------------------------------------------------------------------------
A = 0.15*1e-3; % V
fresp = 0.2667; % 16 resp cycle / min
xyz = [1,0,0.04];
[t,data] = ode45('physiomodel',t,xyz,[],A,fresp,rr,dt,thetai,a,b);
ecg = data(:,3);
save('syn_ecg.mat','ecg','t')