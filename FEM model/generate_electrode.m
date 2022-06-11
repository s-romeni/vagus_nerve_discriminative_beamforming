function model = generate_electrode(model, type, params)
%%-------------------------------------------------------------------------
% Function description:
%%-------------------------------------------------------------------------
% Generate cuff electrode or TIME geometry
%%-------------------------------------------------------------------------
% Inputs:
%%-------------------------------------------------------------------------
% • model: COMSOL object containing the nerve model.
% • type: electrode type.
% • params: electrode geometrical parameters.
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
switch type
    case 'CUFF'
        % Metti l'elettrodo tangente


        R = params(1)/2; % Internal x-length [mm]
        b = params(2); % Internal y-length [mm]
        t = params(3); % Thickness [mm]
        % ell = params(4); % Extrusion length [mm]
        n = params(5); % Number of active sites [1]
        d_as = params(6); % Active site diameter [mm]
        h_as = params(7); % Active site depth [mm]
        if (~isfeasible(params,R))
            error('Immitted data are unfeasible');
        end

        l_cc = 180/(2*n);

        wp1e = geom1.feature.create('wp1e', 'WorkPlane');
        wp1e.set('planetype', 'quick');
        wp1e.set('quickplane', 'xy');
        wp1e.set('quickz',-b/2)

        bodyin = wp1e.geom.feature.create('bodyin','Circle');
        % bodyin.set('base','center');
        bodyin.set('x',0);
        bodyin.set('y',0);
        bodyin.set('r',R);

        geom1.run('wp1e');

        bodyinfull = geom1.feature.create('bodyinfull','Extrude');
        bodyinfull.set('distance',b);

        wp3e = geom1.feature.create('wp3e', 'WorkPlane');
        wp3e.set('planetype', 'quick');
        wp3e.set('quickplane', 'xy');
        wp3e.set('quickz',-b/2)

        bodyout = wp3e.geom.feature.create('bodyout','Circle');
        bodyout.set('x',0);
        bodyout.set('y',0);
        bodyout.set('r',R+t);

        geom1.run('wp3e');

        bodyoutfull = geom1.feature.create('bodyoutfull','Extrude');
        bodyoutfull.set('distance',b);

        geom1.runCurrent

        body = geom1.feature.create('bodyfull', 'Difference');
        body.selection('input').set('bodyoutfull');
        body.selection('input2').set('bodyinfull');
        c = cell(1,n);
        for i = 1:n
            theta = l_cc + (i-1)*180/n;
            O = [R*cos(deg2rad(theta)), R*sin(deg2rad(theta))];
            c{i} = geom1.feature.create(strcat('c',num2str(i)), 'Cylinder');
            c{i}.set('x',O(1));
            c{i}.set('y',O(2));
            c{i}.set('r',d_as/2);
            c{i}.set('h',h_as);
            c{i}.set('axistype', 'spherical');
            c{i}.set('ax2', [90 theta]);
        end
        geom1.runCurrent;
        j = 0;
        cc = cell(1,n);
        s = cell(1,2*n);
        for i = 1:n
            cc{i} = geom1.feature.create(strcat('cc',num2str(i)), 'Mirror');
            cc{i}.selection('input').set(strcat('c',num2str(i)));
            cc{i}.set('keep','on');
            cc{i}.set('pos','0 0 0');
            cc{i}.set('axis','0 1 0');
            j = j+1;
            s{j} = strcat('c',num2str(i));
            j = j+1;
            s{j} = strcat('cc',num2str(i));
        end

        as = geom1.feature.create('as', 'Union');
        as.selection('input').set(s);
        as.set('createselection','on');
        as.set('intbnd','off');

        as_aux = geom1.feature.create('as_aux', 'Copy');
        as_aux.selection('input').set('as');
        as_aux.set('keep','on');

        elec = geom1.feature.create('elec','Difference');
        elec.selection('input').set({'bodyfull'});
        elec.selection('input2').set({'as_aux'});
        elec.set('createselection','on');
        elec.set('intbnd','off');

        geom1.runCurrent
        % mphgeom(model)

case 'TIME'
    l_shaft = params(1);
    w_shaft = params(2);
    h_shaft = params(3);
    l_tip = params(4);
    w_tip = params(5);
    n = params(6);
    d_as = params(7);
    h_as = params(8);

    if (~isfeasible(type, params))
        error('Immitted data are unfeasible');
    end

    l_cc = l_shaft/n;

    %%-----------------------------------------------------------------
    wp1e = geom1.feature.create('wp1e', 'WorkPlane');
    wp1e.set('planetype', 'quick');
    wp1e.set('quickplane', 'xz');

    tab =   [0,         w_shaft/2;  ...
        l_tip,    w_tip/2;    ...
        l_tip,    -w_tip/2;   ...
        0,         -w_shaft/2; ...
        -l_shaft,   -w_shaft/2; ...
        -l_shaft,   w_shaft/2];

    p1 = wp1e.geom.feature.create('p1','Polygon');
    p1.set('type','closed');
    p1.set('source','table');
    p1.set('table',tab);

    p1sol = wp1e.geom.feature.create('p1sol','ConvertToSolid');
    p1sol.selection('input').set({'p1'});

    geom1.run('wp1e');

    shaft = geom1.feature.create('shaft','Extrude');
    shaft.set('distance',-h_shaft/2);
    %%-----------------------------------------------------------------
    wp2e = geom1.feature.create('wp2e', 'WorkPlane');
    wp2e.set('planetype', 'quick');
    wp2e.set('quickplane', 'xz');
    wp2e.set('quicky', h_shaft/2 - h_as);

    for i = 1:n
        O = [-l_cc*(i-1), 0];
        c{i} = wp2e.geom.feature.create(strcat('c',num2str(i)),'Circle');
        c{i}.set('x',O(1));
        c{i}.set('y',O(2));
        c{i}.set('r',d_as/2);
        c{i}.set('type','solid');
    end

    geom1.run('wp2e');

    asup = geom1.feature.create('asup','Extrude');
    asup.set('distance',-h_as);

    as_aux = geom1.feature.create('as_aux','Copy');
    as_aux.selection('input').set('asup');

    elecup = geom1.feature.create('elecup','Difference');
    elecup.selection('input').set({'shaft'});
    elecup.selection('input2').set({'as_aux'});

    geom1.runCurrent
    %%-----------------------------------------------------------------
    elecdwn = geom1.feature.create('elecdwn','Mirror');
    elecdwn.selection('input').set('elecup');
    elecdwn.set('keep','on');
    elecdwn.set('pos','0 0 0');
    elecdwn.set('axis','0 1 0');

    asdwn = geom1.feature.create('asdwn','Mirror');
    asdwn.selection('input').set('asup');
    asdwn.set('keep','on');
    asdwn.set('pos','0 0 0');
    asdwn.set('axis','0 1 0');

    elec = geom1.feature.create('elec','Union');
    elec.selection('input').set({'elecup','elecdwn'});
    elec.set('createselection','on');
    elec.set('intbnd','off');

    as = geom1.feature.create('as','Union');
    as.selection('input').set({'asup','asdwn'});
    as.set('createselection','on');
    as.set('intbnd','off');

    geom1.runCurrent
    % mphgeom(model)

case 'point'
    h = params(1);
    r = params(2);

    as = geom1.feature.create('as','Cylinder');
    as.set('r',r);
    as.set('h',h);
    as.set('createselection','on');

end
end

%%-------------------------------------------------------------------------
% Check the feasibility of the electrode
function bool = isfeasible(type,params)
switch type
    case 'FINE'
        constr1 = (params(3) > params(7));
        constr2 = (params(4) > params(6));
        constr3 = (R > params(5)*params(6));
        bool = constr1 & constr2 & constr3;
    case 'TIME'
        constr1 = (params(2) > params(5));
        constr2 = (params(2) > params(7));
        constr3 = (params(3) > 2*params(8));
        constr4 = (params(1) > params(6)*params(7));
        bool = constr1 & constr2 & constr3 & constr4;
end
end

