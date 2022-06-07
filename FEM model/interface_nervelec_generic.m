function model = interface_nervelec_generic(model, theta, asnum, P, elec_params)
%%-------------------------------------------------------------------------
% Function description: 
%%-------------------------------------------------------------------------
% Interface electrode and nerve geometries
%%-------------------------------------------------------------------------
% Inputs: 
%%-------------------------------------------------------------------------
% • model: COMSOL object containing the nerve model.
% • theta: insertion angle.
% • P: desired location of leading active site.
% • elec_params: electrode geometrical parameters.
%%-------------------------------------------------------------------------
% Outputs:
%%-------------------------------------------------------------------------
% • model: COMSOL object containing the nerve model (updated).
%%-------------------------------------------------------------------------
% Authors: 
%%-------------------------------------------------------------------------
% Andrea Pitzus @TNE, SSSA // @MeDSP, UniCa & Simone Romeni @TNE, EPFL
%%-------------------------------------------------------------------------

geom1 = model.geom('geom1');

% Rename parameters for code readability
l_shaft = elec_params(1);
% w_shaft = elec_params(2);
h_shaft = elec_params(3);
% l_tip = elec_params(4);
% w_tip = elec_params(5);
n = elec_params(6);
% d_as = elec_params (7);
h_as = elec_params(8);
l_cc = l_shaft/n;

% Set leading site
if asnum <= n
    PE = [-l_cc*(asnum-1), h_shaft/2 - h_as/2, 0];
else
    PE = [-l_cc*(asnum - n - 1), -h_shaft/2 + h_as/2, 0];
end
v = P - PE;

% Interface electrode and nerve geometries
rot2 = geom1.feature.create('rot2','Rotate');
rot2.selection('input').set({'elec','as'});
rot2.set('pos', [0, 0, 0]);
rot2.set('axis', [0, 0, 1]);
rot2.set('rot',theta);

mov1 = geom1.feature.create('mov1','Move');
mov1.selection('input').set('rot2');
mov1.set('displ', v);

electrode = geom1.feature.create('electrode','Copy');
electrode.selection('input').set('mov1');

nerve = geom1.feature.create('nerve','Difference');
nerve.selection('input').set({'salgeom','epifull'});
nerve.selection('input2').set('mov1');
end
