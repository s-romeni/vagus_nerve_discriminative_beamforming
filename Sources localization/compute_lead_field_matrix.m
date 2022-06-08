clc
clear
%%-------------------------------------------------------------------------
% General info: Compute the Lead Field Matrix from the FEM model
%%-------------------------------------------------------------------------
% Authors: 
%%-------------------------------------------------------------------------
% Andrea Pitzus @TNE, SSSA // @MeDSP, UniCa & Simone Romeni @TNE, EPFL
%%-------------------------------------------------------------------------
for e = 1:2
    switch e
        case 1
            type = 'TIME';
        case 2
            type = 'CUFF';
    end
    for i_sec = 1:N_sec % current number of nerve model
        %%-----------------------------------------------------------------
        load(['nerve_mod_vagus_human_' num2str(i_sec) '.mat'],'R')
        load(['results_' type '_human_vagus_' num2str(i_sec) '.mat'],'dataref')
        %%-----------------------------------------------------------------
        % Electrodes parameters
        %%-----------------------------------------------------------------
        h_as = 0.03*1e-3;
        r_as = 0.5*0.05*1e-3;
        vol_as = (h_as)*pi*(r_as)^2;
        N_channels = length(dataref);
        %%-----------------------------------------------------------------
        voxels = 40; % 50 micron res
        %%-----------------------------------------------------------------
        x = linspace(-R, R, voxels);
        y = linspace(-R, R, voxels);
        z = 0; % ---> A more accurate lead fiel matrix should take into account
        % at least a portion of the segment not just this value
        [x, y, z] = meshgrid(x, y, z);
        xq = x(:);
        yq = y(:);
        zq = z(:);
        %%-----------------------------------------------------------------
        % Pixels within nerve section
        %%-----------------------------------------------------------------
        V = zeros(round(voxels^2),N_channels);
        for i_channel = 1:N_channels
            F = scatteredInterpolant(dataref{i_channel}.p',dataref{i_channel}.d1');
            vq = F(xq, yq, zq);
            %%-------------------------------------------------------------
            % Matrix assembly
            %%-------------------------------------------------------------
            V(:,i_channel) = vq;
        end
        LFM = (V*(1/vol_as)*eye(N_channels))';
        save(['LFM_' type '_human_vagus_' num2str(i_sec) '.mat'],'LFM')
        %%-----------------------------------------------------------------
    end
end
