clc
clear
%%-------------------------------------------------------------------------
% General info: create a action potential template for each fiber and each
% electrode
%%-------------------------------------------------------------------------
% Along the scripts several parameters have been setted based on the
% following bibliographic source:
% N. Jayaprakash et al., “Organ- and function-specific anatomical
% organization and bioelectronic modulation of the vagus nerve,”
% Neuroscience, preprint, Mar. 2022. doi: 10.1101/2022.03.07.483266.
%%-------------------------------------------------------------------------
% Authors: 
%%-------------------------------------------------------------------------
% Andrea Pitzus @TNE, SSSA // @MeDSP, UniCa & Simone Romeni @TNE, EPFL
%%-------------------------------------------------------------------------
% Analysis are done with a dt = 0.025 ms equally to fs = 40 kHz
%%-------------------------------------------------------------------------
N_sec = 1; % Define the number of section of the model to be generate
%%-------------------------------------------------------------------------
a = [];
for i_sec = 1:N_sec % current number of nerve model
    load(['nerve_mod_vagus_human_' num2str(i_sec) '.mat'])
    N_fasc = size(circular_fascicles,1);
    for i_fasc = 1:N_fasc
        a = [a pi*circular_fascicles(i_fasc, 3)*circular_fascicles(i_fasc, 3)];
    end
end
%%-------------------------------------------------------------------------
meanafffiber = 766;
meanefffiber = 335;
GainAff = 0.25*meanafffiber/(mean(a)); % limit num of fiber to 1/4 to decrease cc
GainEff = 0.25*meanefffiber/(mean(a));
%%-------------------------------------------------------------------------
% Check
%%-------------------------------------------------------------------------
% b = [];
% for i_sec = 1:N_sec % current number of nerve model
%     load(['nerve_mod_vagus_human_' num2str(i_sec) '.mat'])
%     N_fasc = size(circular_fascicles,1);
%     for i_fasc = 1:N_fasc
%         b = [b round(GainAff*pi*circular_fascicles(i_fasc, 3)*circular_fascicles(i_fasc, 3))]; % number of fibers
%     end
% end
%%-------------------------------------------------------------------------

%%-------------------------------------------------------------------------
% Fiber type A-beta (sensory)
%%-------------------------------------------------------------------------
step = 0;
for d = 6:12
    step = step+1;
    load(['Abeta_fiber_' num2str(d, '%1.0d') '.mat'],'data')
    [~,start] = max(data(31,:));
    [~,stop] = max(data(end-30,:));
    data = -1e-3*1e-8*area*data(31:end-30,start-1e2:stop+1e2); % from mA/cm^2 to A
    data = resample(data',1,5)'; % from dt = 0.005 to dt = 0.025 ms
    dx = x(2)-x(1);
    x = [x(1)-dx x x(end)+dx]';
    fiberAb(step).data = data;
    fiberAb(step).x = x(31:end-31);
    fiberAb(step).diameter = d;
    %%---------------------------------------------------------------------
end
%%-------------------------------------------------------------------------
% Fiber type A-alpha (motor)
%%-------------------------------------------------------------------------
step = 0;
for d = 10:13
    step = step+1;
    load(['Aalpha_fiber_' num2str(d, '%1.0d') '.mat'],'data')
    [~,start] = max(data(21,:));
    [~,stop] = max(data(end-20,:));
    data = -1e-3*1e-8*area*data(21:end-20,start-1e2:stop+1e2); % from mA/cm^2 to A
    data = resample(data',1,5)'; % from dt = 0.005 to dt = 0.025 ms
    dx = x(2)-x(1);
    x = [x(1)-dx x x(end)+dx]';
    fiberAa(step).data = data;
    fiberAa(step).x = x(21:end-21);
    fiberAa(step).diameter = d;
    %%---------------------------------------------------------------------
end
%%-------------------------------------------------------------------------
% Pseudo Gaussian Distribution
%%-------------------------------------------------------------------------
dst_beta = round(sqrt(8)*randn(1e6,1)+4); % mean diameter = 9 um
dst_beta(dst_beta < 1 | dst_beta > 7) = [];
%%-------------------------------------------------------------------------
dst_alpha = round(sqrt(5)*randn(1e6,1)+2.5); % mean diameter = 11.5 um
dst_alpha(dst_alpha < 1 | dst_alpha > 4) = [];
%%-------------------------------------------------------------------------
dt = 0.025*1e-3; % time resolution
%%-------------------------------------------------------------------------
for e = 1:2
    switch e
        case 1
            type = 'TIME';
        case 2
            type = 'CUFF';
    end
    h_as = 0.3*1e-3;
    r_as = 0.5*0.5*1e-3;
    vol_as = (h_as)*pi*(r_as)^2;
    %%---------------------------------------------------------------------
    for i_sec = 1:N_sec % current number of nerve model
        load(['nerve_mod_vagus_human_' num2str(i_sec) '.mat'],'circular_fascicles','circular_fascicles_TIME','R')
        load(['results_' type '_human_vagus_' num2str(i_sec) '.mat'],'dataref')
        %%-----------------------------------------------------------------
        circular_fascicles_CUFF = circular_fascicles;
        if strcmp(type,'TIME')
            circular_fascicles = circular_fascicles_TIME;
        end
        N_channels = length(dataref);
        %%-----------------------------------------------------------------
        % Nerve topography
        %%-----------------------------------------------------------------
        if strcmp(type,'TIME')
            [~,I] = sort(circular_fascicles(:,3));
            candidate_motor = I(1:end-1);
            motor_fasc = candidate_motor(randi([1 length(candidate_motor)],1));
            defferent = zeros(fasc_number,1);
            for i_fib = 1:fasc_number
                defferent(i_fib) = norm(circular_fascicles_CUFF(motor_fasc,[1 2]) - circular_fascicles_CUFF(i_fib,[1 2]));
            end
            [~,I] = sort(defferent);
            motor_fasc = I(1:floor(0.5*fasc_number));
        else
            load(['nerve_mod_vagus_human_' num2str(i_sec) '_TIME_rec.mat'],'motor_fasc')
        end
        %%-----------------------------------------------------------------
        % Assign fiber location and size
        %%-----------------------------------------------------------------
        N_fasc = size(circular_fascicles,1);
        voxels = 2000; % 1 micron resolution on human vagus
        xg = linspace(-R, R, voxels);
        yg = linspace(-R, R, voxels);
        z = 0;
        [x, y, zq] = meshgrid(xg, yg, z);
        fibers = struct('location',[],'data',[]);
        fascicles = cell(N_fasc,1);
        for i_fasc = 1:N_fasc
            sel = (x-circular_fascicles(i_fasc, 1)).^2 + (y-circular_fascicles(i_fasc, 2)).^2 < circular_fascicles(i_fasc, 3)^2;
            [row,col] = find(sel);
            if sum(fasc_index == motor_fasc) > 0
                for i_fib = 1:round(GainEff*circular_fascicles(fasc_index, 3)*circular_fascicles(fasc_index, 3)) % number of fibers
                    % Num of fibers = GainEff*[0.15, 0.7]*1e-3).^2
                    % Area fascicle = (([0.15, 0.7]*1e-3).^2)*pi
                    nrand = randi(length(row));
                    pos = [x(row(nrand),col(nrand)),y(row(nrand),col(nrand))];
                    diam = datasample(dst_alpha,1,'Replace',false);
                    diameter = fiberAa(diam).diameter;
                    sel = sel & ((x-pos(1)).^2 + (y-pos(2)).^2 > (0.5*1e-6*diameter)^2);
                    [row,col] = find(sel);
                    if length(row) >= 1
                        nrand = randi(length(row));
                        fibers(i_fib).location(1:2) = pos;
                        fibers(i_fib).location(3) = diameter;
                        Re = zeros(N_channels,length(fiberAa(diam).x));
                        for i_channel = 1:N_channels
                            %%---------------------------------------------
                            % FEM Analysis
                            %%---------------------------------------------
                            F = scatteredInterpolant(dataref{i_channel}.p',dataref{i_channel}.d1');
                            [~, ~, zg] = meshgrid(fibers(i_fib).location(1), fibers(i_fib).location(2), fiberAa(diam).x);
                            zq = zg(:);
                            r = zeros(1,length(fiberAa(diam).x));
                            for node = 1:length(zq)
                                vq = F(fibers(i_fib).location(1), fibers(i_fib).location(2), zq(node));
                                Res = vq/vol_as;
                                r(1,node) = Res;
                            end
                            Re(i_channel,:) = r;
                        end
                        % Re = fliplr(Re); % efferent (since is symmetric
                        % there is no need)
                        fibers(i_fib).data = Re*fiberAa(diam).data;
                    end
                end
            else
                for i_fib = 1:round(GainAff*circular_fascicles(fasc_index, 3)*circular_fascicles(fasc_index, 3))
                    nrand = randi(length(row));
                    pos = [x(row(nrand),col(nrand)),y(row(nrand),col(nrand))];
                    diam = datasample(dst_beta,1,'Replace',false);
                    diameter = fiberAb(diam).diameter;
                    sel = sel & ((x-pos(1)).^2 + (y-pos(2)).^2 > (0.5*1e-6*diameter)^2);
                    [row,col] = find(sel);
                    if length(row) >= 1
                        nrand = randi(length(row));
                        fibers(i_fib).location(1:2) = pos;
                        fibers(i_fib).location(3) = diameter;
                        Re = zeros(N_channels,length(fiberAb(diam).x));
                        for i_channel = 1:N_channels
                            %%---------------------------------------------
                            % FEM Analysis
                            %%---------------------------------------------
                            F = scatteredInterpolant(dataref{i_channel}.p',dataref{i_channel}.d1');
                            [~, ~, zg] = meshgrid(fibers(i_fib).location(1), fibers(i_fib).location(2), fiberAb(diam).x);
                            zq = zg(:);
                            r = zeros(1,length(fiberAb(diam).x));
                            for node = 1:length(zq)
                                vq = F(fibers(i_fib).location(1), fibers(i_fib).location(2), zq(node));
                                Res = vq/vol_as;
                                r(1,node) = Res;
                            end
                            Re(i_channel,:) = r;
                        end
                        % Re = Re; % afferent
                        fibers(i_fib).data = Re*fiberAb(diam).data;
                    end
                end
            end
            fascicles{i_fasc} = fibers;
            clear fibers
        end
        nerve.fascicles = fascicles;
        save(['nerve_mod_vagus_human_' num2str(i_sec) '_' type '_rec.mat'],'nerve','motor_fasc')
    end
end
