# Discriminative Beamforming
This is the repository for the manuscript *A method to establish functional vagus nerve topography from electro-neurographic spontaneous activity*, under review at [Cell: Patterns](https://www.cell.com/patterns/home). The manuscript pre-print is available [here].

**Authors:** Andrea Pitzus*, Simone Romeni*, Fabio Vallone, Silvestro Micera.
**equally contribute*
## A short guide to the repositories

![Workflow](readme.png)

### Main scripts
These repo contains the following scripts:
##### [FEM model](https://github.com/s-romeni/vagus_nerve_discriminative_beamforming/tree/main/FEM%20model) 
* `fem_model_vagus_nerve.m`, which generate a plausible human left vagus nerve model.
##### [Record action potentials](https://github.com/s-romeni/vagus_nerve_discriminative_beamforming/tree/main/Record%20action%20potentials) 
* `Aalpha_fibers_nodes_current.py`,`Abeta_fibers_nodes_current.py`, which generate the nodes currents of myelinated motor and sensory fibers;
* `action_potential_templates.m`, which simulates the action potential recording of each fiber.
##### [Generate electro-neurographic data](https://github.com/s-romeni/vagus_nerve_discriminative_beamforming/tree/main/Generate%20electro-neurographic%20data) 
* `generate_blood_pressure.m`,`generate_pulmonary_volume.m`, which generate the simulated average blood pressure and pulmonary volume in resting conditions;
* `generate_eng_data.m`, which generate the simulated electro-neurographic signals recorded from the human left vagus nerve.
##### [Sources localization](https://github.com/s-romeni/vagus_nerve_discriminative_beamforming/tree/main/Sources%20localization) 
* `compute_lead_field_matrix.m`, which compute the lead field matrix from the human left vagus nerve model;
* `single_source_localization.m`, which perform the single source localization case;
* `multiple_sources_localization.m`, which perform the multiple sources localization case;
* `biophysical_sources_localization.m`, which perform the biophysical sources localization case-

All data needed to execute the scripts are contained in the folder.
