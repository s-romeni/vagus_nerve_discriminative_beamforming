# Establish vagus nerve functional topography using different source localization methods
This is the repository for the manuscript *A method to establish functional vagus nerve topography from electro-neurographic spontaneous activity*, under review at [Cell: Patterns](https://www.cell.com/patterns/home). The manuscript pre-print is available [here].

**Authors:** Andrea Pitzus*, Simone Romeni*, Fabio Vallone, Silvestro Micera.
**equally contribute*
## A short guide to the repositories

![Workflow](readme.png)

### Main scripts
##### [:file_folder: FEM model](https://github.com/s-romeni/vagus_nerve_discriminative_beamforming/tree/main/FEM%20model) 
* `fem_model_vagus_nerve.m`, which generate a plausible human left vagus nerve model.
##### [:file_folder: Record action potentials](https://github.com/s-romeni/vagus_nerve_discriminative_beamforming/tree/main/Record%20action%20potentials) 
* `Aalpha_fibers_nodes_current.py`,`Abeta_fibers_nodes_current.py`, which generate the nodes currents of myelinated motor and sensory fibers;
* `action_potential_templates.m`, which simulates the action potential recording of each fiber.
##### [:file_folder: Generate electro-neurographic data](https://github.com/s-romeni/vagus_nerve_discriminative_beamforming/tree/main/Generate%20electro-neurographic%20data) 
* `generate_blood_pressure.m`,`generate_pulmonary_volume.m`, which generate the simulated average blood pressure and pulmonary volume in resting conditions;
* `generate_eng_data.m`, which generate the simulated electro-neurographic signals recorded from the human left vagus nerve.
##### [:file_folder: Sources localization](https://github.com/s-romeni/vagus_nerve_discriminative_beamforming/tree/main/Sources%20localization) 
* `compute_lead_field_matrix.m`, which compute the lead field matrix from the human left vagus nerve model;
* `single_source_localization.m`, which perform the single source localization case;
* `multiple_sources_localization.m`, which perform the multiple sources localization case;
* `biophysical_sources_localization.m`, which perform the biophysical sources localization case.

### Other functions
##### [`fem_model_vagus_nerve.m`]
* `[centers,r] = aci_packing(R, rmax, rmin, delta, epsilon, N)`, A-priori Check for Intersections packing of circular objects in a circular region.
* `circular_fascicles = reshape_nerve(circular_fascicles,h,R,delta)`, iteratevely move the fascicles to mimic the nerve re-organization after the implantation of a intraneural interface, in particular the TIME.  
* `draw_section(R,circular_fascicles_TIME,circular_fascicles,l_shaft,pos,d_as,h_as)`, draw the nerve section to check the topography and the reorganization. 
* `model = generate_ec(model)`, sets up the physical interface ec (as a Conductive Media) and saves the type of solver that will be used.
* `model = generate_outernerve(model, type, opt, ell, epineurium, sal_pars, custom)`, generate saline bath and epineurium. 
* `model = generate_circfasc(model, circular_fascicles, ell)`, generate elliptical fascicles geometry. 
* `model = generate_electrode(model, type, params)`, generate cuff electrode or TIME geometry. 
* `model = interface_nervelec(model, theta, asnum, P, elec_params)`, interface electrode and nerve geometries. 
* `model = interface_nervelec_generic(model, theta, asnum, P, elec_params)`, interface electrode and nerve geometries in the generic model case. 
* `assign_materials(model)`, assign materials to the different domains. 
* `assign_materials_generic(model)`, assign materials to the different domains in the generic model case. 
##### [`generate_blood_pressure.m`]
* `data = physiomodel(t,xyz,A,fresp,rr,dt,thetai,a,b)`, generate synthetic physiological data given the model parameters.



All data needed to execute the scripts are generated within the script.
