clc
clear
close all
%%-------------------------------------------------------------------------
% General info: create a FEM model of a plausible human vagus nerve
%%-------------------------------------------------------------------------
% Along the scripts several parameters have been setted based on the  
% following bibliographic source:
% N. A. Pelot et al., "Quantified Morphology of the Cervical and 
% Subdiaphragmatic Vagus Nerves of Human, Pig, and Rat,” Front. Neurosci., 
% vol. 14, p. 601479, Nov. 2020, doi: 10.3389/fnins.2020.601479.
%%-------------------------------------------------------------------------
import com.comsol.model.*
import com.comsol.model.util.*
%%-------------------------------------------------------------------------
N_sec = 1; % Define the number of section of the model to be generate
nerve_diam = 2; % [mm]
for i_sec = 1:N_sec % current number of nerve model 
    %%---------------------------------------------------------------------
    % Set TIME (electrode) parameters
    %%---------------------------------------------------------------------
    % elec_params = [Shaft length [mm], Shaft width [mm], Shaft height [mm], ...
    % Tip length [mm],Tip width [mm], Num of channels x side [1], ... 
    % Active site diameter [mm], Active site height [mm], ...
    % Inter channels distance [mm]];
    elec_params = [nerve_diam-(0.20+0.05), 0.5, 0.1, 0.20, 0.3, 7, 0.05, 0.03];
    elec_params([1 2 3 4 5 7 8]) = elec_params([1 2 3 4 5 7 8]).*1e-3; % mm2m conversion
    l_shaft = elec_params(1);
    w_shaft = elec_params(2);
    h_shaft = elec_params(3);
    l_tip   = elec_params(4);
    w_tip   = elec_params(5);
    N_asxs  = elec_params(6);
    d_as    = elec_params(7);
    h_as    = elec_params(8);
    l_cc    = l_shaft/N_asxs;
    %%---------------------------------------------------------------------
    % Set nerve section parameters
    %%---------------------------------------------------------------------
    % nerve_pars ['Number of fascicles [1]', 'Fascicle minimum radius [mm]', ...
    %     'Fascicle maximum radius [mm]', 'Nerve external radius [mm]', ...
    %     'delta [mm]', 'epsilon [mm]', 'Nerve extrusion length [mm]'};
    nerve_pars = [6, 0.08, 0.35, nerve_diam/2, 0.05, 0.01, 20];
    nerve_pars(2:end) = nerve_pars(2:end)*1e-3; % mm2m conversion
    N_fasc  = nerve_pars(1);
    rmin    = nerve_pars(2);
    rmax    = nerve_pars(3);
    R       = nerve_pars(4);
    delta   = nerve_pars(5);
    epsilon = nerve_pars(6);
    %%---------------------------------------------------------------------
    % Set saline parameters
    sal_pars = [4*1e-3, 6*1e-3]; % [Radius [mm], + Delta Lenght [mm]]    
    %%---------------------------------------------------------------------
    [centers,radii] = aci_packing(R, rmax, rmin, delta, epsilon, N_fasc); % A-priori Check for Intersections (no overlapped fascicles)
    circular_fascicles = [centers, radii];
    %%---------------------------------------------------------------------
    % Modelling nerve reorganization
    %%---------------------------------------------------------------------
    % Move fascicles above the electrode
    circular_fascicles_TIME = reshape_nerve(circular_fascicles,(-h_shaft/2-h_as+h_shaft)+epsilon,R,epsilon);
    % Move fascicles below the electrode
    circular_fascicles_TIME = reshape_nerve(circular_fascicles_TIME,-h_shaft/2-h_as-epsilon,R,epsilon);
    %%---------------------------------------------------------------------
    % Active site position
    %%---------------------------------------------------------------------
    r_as = d_as/2;
    pos = zeros(2*N_asxs,2);
    for i = 1:2*N_asxs
        if i <= N_asxs
            pos(i,:) = [-l_cc*(i-1) + r_as + (nerve_diam/2 - (tl+0.05))*1e-3, h_shaft/2];
        else
            pos(i,:) = [-l_cc*(i-N_asxs-1) + r_as + (nerve_diam/2 - (tl+0.05))*1e-3, -h_shaft/2];
        end
    end
    %%---------------------------------------------------------------------
    % Draw the nerve section to check the topography and the reorganization
    %%---------------------------------------------------------------------
    draw_section(R,circular_fascicles_TIME,circular_fascicles,l_shaft,pos,d_as,h_as);
    %%---------------------------------------------------------------------
    if not(isfile(['nerve_mod_vagus_human_' num2str(i_sec) '.mat']))
        %%-----------------------------------------------------------------
        save(['nerve_mod_vagus_human_' num2str(i_sec) '.mat'], 'R', 'ell','circular_fascicles','circular_fascicles_TIME');
    else
        %%-----------------------------------------------------------------
        load(['nerve_mod_vagus_human_' num2str(i_sec) '.mat'], 'R', 'ell','circular_fascicles','circular_fascicles_TIME');
    end
    %%---------------------------------------------------------------------
    % Create model and assign conditions
    %%---------------------------------------------------------------------
    dataref = cell(2*N_asxs,1);
    channel = 0;
    h = waitbar(channel/2*N_asxs,'Please wait...');
    for channel = 1:2*N_asxs
        model = ModelUtil.create('Model');
        ModelUtil.showProgress(true); % No 'component' here
        geom1 = model.geom.create('geom1', 3); % 3 stay for 3D (LiveLink Guide pag 58)
        model = generate_ec(model); % Create physics
        %%-----------------------------------------------------------------

        %%-----------------------------------------------------------------
        % Interface nerve and electrode models
        %%-----------------------------------------------------------------
        model = generate_outernerve(model, 'TIME', [0, 0, 0], ell, R, sal_pars, 0);
        %%-----------------------------------------------------------------
        model = generate_circfasc(model, circular_fascicles_TIME, ell);
        %%-----------------------------------------------------------------
        model = generate_electrode(model, 'TIME', elec_params);
        %%-----------------------------------------------------------------
        % Select boundaries of active site
        %%-----------------------------------------------------------------
        delta = 1e-6;
        idx = cell(2*N_asxs,1);
        for i = 1:2*N_asxs
            if i <= N_asxs
                p0 = [-l_cc*(i-1) + d_as/2 + delta, ...
                    h_shaft/2 - h_as - delta,         ...
                    -d_as/2 - delta];
                p1 = [-l_cc*(i-1) - d_as/2 - delta, ...
                    h_shaft/2 + delta,                ...
                    d_as/2 + delta];
            else
                p0 = [-l_cc*(i - N_asxs - 1) + d_as/2 + delta,   ...
                    -h_shaft/2 + h_as + delta,                ...
                    -d_as/2 - delta];
                p1 = [-l_cc*(i - N_asxs - 1) - d_as/2 - delta,   ...
                    -h_shaft/2 - delta,                       ...
                    d_as/2 + delta];
            end
            idx{i} = mphselectbox(model,'geom1',[p0',p1'],'domain');
        end
        %%-----------------------------------------------------------------
        curr = zeros(2*N_asxs,1);
        curr(channel) = 1; % [A/m^3]
        %%-----------------------------------------------------------------
        for i = 1:2*N_asxs
            model.physics('ec').feature.create(strcat('cs', num2str(i)), 'CurrentSource', 3);
            model.physics('ec').feature(strcat('cs', num2str(i))).selection.set(idx{i});
            model.physics('ec').feature(strcat('cs', num2str(i))).set('Qj', curr(i));
        end
        %%-----------------------------------------------------------------
        % Set insertion parameters
        %%-----------------------------------------------------------------
        % insertion_params = [Leading active site [1], Desired x-position of active site [mm]', ...
        % 'Desired y-position of active site [mm]', 'Desired angle of insertion [deg]'];
        insertion_params = [1 l_shaft 0 0];
        asnum = insertion_params(1);
        P(1)  = insertion_params(2);
        P(2)  = insertion_params(3);
        P(3)  = 0;
        theta = insertion_params(4);
        %%-----------------------------------------------------------------
        model = interface_nervelec(model, theta, asnum, P, elec_params);
        model.geom('geom1').run;
        %%-----------------------------------------------------------------
        % Material
        %%-----------------------------------------------------------------
        model = generate_materials(model);
        model = assign_materials(model);
        %%-----------------------------------------------------------------
        % Geometry
        %%-----------------------------------------------------------------
        model.geom('geom1').run;
        %%-----------------------------------------------------------------
        % Mesh
        %%-----------------------------------------------------------------
        model.mesh.create('mesh1', 'geom1');
        model.mesh('mesh1').autoMeshSize(3);
        model.mesh('mesh1').run;
        stats = mphmeshstats(model);
        ModelUtil.showProgress(true)
        %%-----------------------------------------------------------------
        % Run
        %%-----------------------------------------------------------------
        model.study('std1').run;
        dataref{channel} = mpheval(model,{'V'});
        waitbar(channel/(2*N_asxs),h,[num2str(channel) '/' num2str(2*N_asxs)])
        %%-----------------------------------------------------------------
    end
    close(h)
    save(['results_TIME_human_vagus_' num2str(i_sec) '.mat'],'dataref')
    %%---------------------------------------------------------------------
    
    %%---------------------------------------------------------------------
    % Set CUFF electrode parameters
    %%---------------------------------------------------------------------
    % elec_params = [Internal length [mm], Internal z-length [mm], Thickness [mm], ...
    % Extrusion length [mm], Num of channels x side [1], ... 
    % Active site diameter [mm], Active site height [mm], Active site diameter [mm], ...
    % Active site depth [mm]];
  
    %%---------------------------------------------------------------------
    elec_params = [nerve_diam, 0.5, 0.1, 10, 7, 0.05, 0.03];
    elec_params([1 2 3 4 6 7]) = elec_params([1 2 3 4 6 7]).*1e-3; % mm2m
    a    = elec_params(1);
    b    = elec_params(2);
    t    = elec_params(3);
    % ell  = elec_params(4);
    N_asxs = elec_params(5);
    d_as = elec_params(6);
    h_as = elec_params(7);
    l_cc = 180/(2*N_asxs);
    opt = [0, elec_params(1), elec_params(2)];
    %%---------------------------------------------------------------------
    dataref = cell(2*N_asxs,1);
    channel = 0;
    h = waitbar(channel/2*N_asxs,'Please wait...');
    for channel = 1:2*N_asxs
        model = ModelUtil.create('Model');
        ModelUtil.showProgress(true); % No 'component' here
        geom1 = model.geom.create('geom1', 3); % 3 stay for 3D (LiveLink Guide pag 58)
        model = generate_ec(model); % Create physics
        %%-----------------------------------------------------------------

        %%-----------------------------------------------------------------
        model = generate_outernerve(model, 'CUFF', opt, ell, R, sal_pars, 0);
        %%-----------------------------------------------------------------
        model = generate_circfasc(model, circular_fascicles, ell);
        %%-----------------------------------------------------------------
        model = generate_electrode(model, 'CUFF', elec_params);
        %%-----------------------------------------------------------------
        % Select boundaries of active site
        %%-----------------------------------------------------------------
        delta = 3e-4;
        O = zeros(2*N_asxs,6);
        idx = zeros(2*N_asxs,1);
        for i = 1:2*N_asxs
            if i <= N_asxs
                theta = l_cc + (i-1)*180/N_asxs;
                O(i,1:3) = [R*cos(deg2rad(theta))-delta, R*sin(deg2rad(theta))-delta/1.5,-delta];
                O(i,4:6) = [R*cos(deg2rad(theta))+delta, R*sin(deg2rad(theta))+delta/1.5,+delta];
            else
                theta = pi + l_cc + (i-1)*180/N_asxs;
                O(i,1:3) = [R*cos(deg2rad(theta))-delta, R*sin(deg2rad(theta))-delta,-delta];
                O(i,4:6) = [R*cos(deg2rad(theta))+delta, R*sin(deg2rad(theta))+delta,+delta];
            end
            idx(i) = mphselectbox(model,'geom1',[O(i,1:3)',O(i,4:6)'],'domain');
        end
        %%-----------------------------------------------------------------
        curr = zeros(2*N_asxs,1);
        curr(channel) = 1;
        %%-----------------------------------------------------------------
        for i = 1:2*N_asxs
            model.physics('ec').feature.create(strcat('cs', num2str(i)), 'CurrentSource', 3);
            model.physics('ec').feature(strcat('cs', num2str(i))).selection.set(idx(i));
            model.physics('ec').feature(strcat('cs', num2str(i))).set('Qj', curr(i));
        end
        %%-----------------------------------------------------------------
        % Material
        %%-----------------------------------------------------------------
        model = generate_materials(model);
        model = assign_materials(model);
        %%-----------------------------------------------------------------
        % Interface nerve and electrode models
        %%-----------------------------------------------------------------
        eleccopy = geom1.feature.create('eleccopy', 'Copy');
        eleccopy.selection('input').set('elec');
        %%-----------------------------------------------------------------
        ascopy = geom1.feature.create('ascopy', 'Copy');
        ascopy.selection('input').set('as');
        %%-----------------------------------------------------------------
        nerve = geom1.feature.create('nerve','Difference');
        %%-----------------------------------------------------------------
        nerve.selection('input').set({'salgeom','epigeom','perigeom','endogeom'});
        %%-----------------------------------------------------------------
        nerve.selection('input2').set({'eleccopy','ascopy'});
        %%-----------------------------------------------------------------
        % Geometry
        %%-----------------------------------------------------------------
        % model.geom('geom1').run;
        %%-----------------------------------------------------------------
        % Mesh
        %%-----------------------------------------------------------------
        model.mesh.create('mesh1', 'geom1');
        model.mesh('mesh1').autoMeshSize(3);
        model.mesh('mesh1').run;
        stats = mphmeshstats(model);
        ModelUtil.showProgress(true)
        %%-----------------------------------------------------------------
        % Run
        %%-----------------------------------------------------------------
        model.study('std1').run;
        dataref{channel} = mpheval(model,{'V'});
        waitbar(channel/(2*N_asxs),h,[num2str(channel) '/' num2str(2*N_asxs)])
    end
    close(h)
    save(['results_CUFF_human_vagus_' num2str(i_sec) '.mat'],'dataref')
end
