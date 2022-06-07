function model = assign_materials_apriori(model)
%%-------------------------------------------------------------------------
% Function description: 
%%-------------------------------------------------------------------------
% Assign materials to the different domains
%%-------------------------------------------------------------------------
% Inputs: 
%%-------------------------------------------------------------------------
% • model: COMSOL object containing the nerve model.
%%-------------------------------------------------------------------------
% Outputs:
%%-------------------------------------------------------------------------
% • model: COMSOL object containing the nerve model (updated).
%%-------------------------------------------------------------------------
% Authors: 
%%-------------------------------------------------------------------------
% Andrea Pitzus @TNE, SSSA // @MeDSP, UniCa & Simone Romeni @TNE, EPFL
%%-------------------------------------------------------------------------

salsel = mphgetselection(model.selection('geom1_salgeom_dom'));
saldom = salsel.entities;
model.material('sal').selection.set(saldom);

episel = mphgetselection(model.selection('geom1_epifull_dom'));
% episel = mphgetselection(model.selection('geom1_epigeom_dom'));
epidom = episel.entities;
model.material('epi').selection.set(epidom);

% perisel = mphgetselection(model.selection('geom1_perigeom_dom'));
% peridom = perisel.entities;
% model.material('peri').selection.set(peridom);
% 
% endosel = mphgetselection(model.selection('geom1_endogeom_dom'));
% endodom = endosel.entities;
% model.material('endo').selection.set(endodom)

elecsel = mphgetselection(model.selection('geom1_elec_dom'));
elecdom = elecsel.entities;
model.material('elec').selection.set(elecdom)

assel = mphgetselection(model.selection('geom1_as_dom'));
asdom = assel.entities;
model.material('as').selection.set(asdom)
end