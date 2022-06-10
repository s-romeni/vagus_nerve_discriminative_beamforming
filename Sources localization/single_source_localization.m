clc
clear
%%-------------------------------------------------------------------------
% General info: Perform a single source localization
%%-------------------------------------------------------------------------
% Along the scripts several methods have been implemented/developed based 
% on the following bibliographic sources:
% B. Wodlinger and D. M. Durand, “Localization and Recovery of Peripheral 
% Neural Sources With Beamforming Algorithms,” IEEE Transactions on Neural 
% Systems and Rehabilitation Engineering, vol. 17, no. 5, pp. 461–468, 
% Oct. 2009, doi: 10.1109/TNSRE.2009.2034072.
% F. Babiloni et al., “Multimodal integration of high-resolution EEG and
% functional magnetic resonance imaging data: a simulation study,” 
% NeuroImage, vol. 19, no. 1, pp. 1–15, May 2003, 
% doi: 10.1016/S1053-8119(03)00052-1.
%%-------------------------------------------------------------------------
% Authors: 
%%-------------------------------------------------------------------------
% Andrea Pitzus @TNE, SSSA // @MeDSP, UniCa & Simone Romeni @TNE, EPFL
%%-------------------------------------------------------------------------
% Analysis are done with a dt = 0.01 s equally to fs = 100 Hz
%%-------------------------------------------------------------------------
% Challenges:
% Is performed a challenge for each fascicle (at time)
% For each challenge the source signal is imposed to be a sin
% All fascicles are a sorce of brownian noise
%%-------------------------------------------------------------------------
N_sec = 1; % Define the number of section of the model to be generate
Model = 'Structural'; % 'Generic' otherwise
%%-------------------------------------------------------------------------

%%-------------------------------------------------------------------------
for e = 1:2
    switch e
        case 1
            type = 'TIME';
            s = 1;
        case 2
            type = 'CUFF';
            s = 1;
    end
    %%---------------------------------------------------------------------
    load('nerve_mod_vagus_human_generic.mat')
    voxels = 40; % 50 micron res
    x = linspace(-R, R, voxels);
    y = linspace(-R, R, voxels);
    z = 0;
    [x, y, z] = meshgrid(x, y, z);
    sel_gen = x.^2 + y.^2 < R^2;
    if strcmp('TIME',type)
        tl = 0.2;
        nerve_diam = 2*R*1e3; % human
        l_shaft = (nerve_diam-(tl+0.05))*1e-3;
        % w_shaft = 1e-3;
        h_shaft = 1e-4;
        h_as_t = 0.03*1e-3;
        r_as_t = 0.5*0.05*1e-3;
        % vol_as = (h_as)*pi*(r_as)^2;
        Nasxs = 7;
        l_cc = l_shaft/Nasxs;
        sel_gen = sel_gen - (x > (-1e-3*nerve_diam/2+l_cc/2) & x < (-1e-3*nerve_diam/2+l_cc/2+l_shaft) & y > (-h_shaft/2-h_as_t) & y < (-h_shaft/2-h_as_t+h_shaft));
    end
    %---------------------------------------------------------------------
    
    %%---------------------------------------------------------------------
    for i_sec = 1:N_sec % current number of nerve model
        %%-----------------------------------------------------------------
        % Integrate structural information
        %%-----------------------------------------------------------------
        load(['nerve_mod_vagus_human_' num2str(i_sec) '.mat'],'circular_fascicles','circular_fascicles_TIME')
        N_fasc = size(circular_fascicles,1);
        if strcmp(type,'TIME')
            circular_fascicles = circular_fascicles_TIME;
        end
        sel = sel_gen;
        if strcmp(model,'Structural')
            for i_fasc = 1:N_fasc
                sell = (x-circular_fascicles(i_fasc, 1)).^2 + (y-circular_fascicles(i_fasc, 2)).^2 <(circular_fascicles(i_fasc, 3)).^2;
                sel = sel+sell;
            end 
        end
        N_vox = voxels^2;
        sel = reshape(sel,[1 N_vox]);
        Winv = diag(sel(:));
        if strcmp(model,'Structural')
            load(['LFM_' type '_human_vagus_' num2str(i_sec) '.mat'],'LFM')
        else
            load(['LFM_' type '_human_vagus_generic.mat'],'LFM')
        end
        hat_LFM = (Winv*LFM')';
        %%-----------------------------------------------------------------
        % Compute the pseudoinverse of the lead field matrix
        %%-----------------------------------------------------------------
        LFM_psd_inv = (LFM*Winv*LFM')^-1;
        LFM_psd_inv = LFM_psd_inv*LFM*Winv;
        LFM_psd_inv = LFM_psd_inv';
        %%-----------------------------------------------------------------
        for i_fasc = 1:size(LFM_psd_inv,1)
            if norm(Winv*LFM'*LFM_psd_inv(i_fasc,:)') > 0
                LFM_psd_inv(i_fasc,:) = LFM_psd_inv(i_fasc,:)/norm(Winv*LFM'*LFM_psd_inv(i_fasc,:)');
            end
        end
        %%-----------------------------------------------------------------

        %%-----------------------------------------------------------------
        N_ch = size(LFM_psd_inv,2);
        N_vox = size(LFM_psd_inv,1);
        %%-----------------------------------------------------------------

        %%-----------------------------------------------------------------
        dt = 0.01;
        fs = 1/dt;
        t = (dt:dt:300);
        %%-----------------------------------------------------------------
        A = 1e-6; % 1 mA
        f = 4; % Hz
        sig = A*sin(2*pi*f*t);
        for snr = 1
            %%-------------------------------------------------------------
            % Generate the signals
            %%-------------------------------------------------------------
            rest_sig = zeros(N_vox,length(t));
            SNR = zeros(N_fasc,N_ch);
            rec_sources = cell(N_fasc,1);
            rec_noise = cell(N_fasc,1);
            for i_fasc = 1:N_fasc
                %%---------------------------------------------------------
                % All the fascicle
                pixs = (x-circular_fascicles(i_fasc, 1)).^2 + (y-circular_fascicles(i_fasc, 2)).^2 < circular_fascicles(i_fasc, 3)^2;
                pixs = reshape(pixs,[1 N_vox]);
                pixs = find(pixs == 1);
                %%---------------------------------------------------------
                rng(s+i_fasc,'twister');
                cn = dsp.ColoredNoise('brown','SamplesPerFrame',length(sig),'NumChannels',length(sel)); % Flicker noise
                rng(s+i_fasc,'twister');
                rest_sig(pixs,:) = cn()';
                %%---------------------------------------------------------
                % Only the central pixel
                dist = (x-circular_fascicles(i_fasc, 1)).^2 + (y-circular_fascicles(i_fasc, 2)).^2;
                dist = reshape(dist,[1 N_vox]);
                [h,I] = sort(dist,'ascend');
                pix = I(1);
                %%---------------------------------------------------------
                rest_sig(pix,:) = zeros(length(pix),1)*sig;
                rec_rest = LFM*(rest_sig);
                S = zeros(size(rest_sig));
                S(pix,:) = ones(length(pix),1)*sig;
                rec_sig = LFM*S;
                rec_sources{i_fasc} = rec_sig;
                rec_noise{i_fasc} = rec_rest;
                %%---------------------------------------------------------
                % Signal-to-Noise Ratio
                %%---------------------------------------------------------
                SNR(i_fasc,:) = rms(rec_sig')./rms(rec_rest');
                %%---------------------------------------------------------
            end
            %%-------------------------------------------------------------
            % mod_noise A = noise amplitude
            mod_noise = mean(mean(SNR));
            %%-------------------------------------------------------------
            % Recordings
            %%-------------------------------------------------------------
            N_challenges = N_fasc;
            recording = cell(N_challenges,1);
            estimated_sources = cell(N_challenges,1);
            SNR = zeros(N_challenges,N_ch);
            for i_challenge = 1:N_challenges
                recording{i_challenge} = rec_sources{i_challenge}+(1/snr)*mod_noise*rec_noise{i_challenge};
                %%---------------------------------------------------------
                % Signal-to-Noise Ratio
                %%---------------------------------------------------------
                SNR(i_challenge,:) = rms(rec_sources{i_challenge}')./rms((1/snr)*mod_noise*rec_noise{i_challenge}');
                %%---------------------------------------------------------
                % Source localization
                estimated_sources{i_challenge} = LFM_psd_inv*recording{i_challenge};
            end
            check = mean(mean(SNR)); % should be = 1
            mod_noise = (1/snr)*mod_noise;
            save(['SNR = ' num2str(snr) '\mod_noise_' type '_' num2str(i_sec) '.mat'],'mod_noise')
            %%-------------------------------------------------------------
            
            %%-------------------------------------------------------------
            % Discriminative Indexes
            %%-------------------------------------------------------------
            N = length(t);
            frequency = ((1:N)-N/2)*fs/N;
            D = zeros(N_challenges,N_ch);
            for i_challenge = 1:N_challenges
                challenge = recording{i_challenge};
                for i_ch = 1:N_ch
                    spectrum = abs(fftshift(fft(challenge(i_ch,:))))/N;
                    D(i_challenge,i_ch) = 2*sum(spectrum(frequency > f-0.1 & frequency < f+0.1));
                end
            end
            %%-------------------------------------------------------------

            %%-------------------------------------------------------------
            % Beamfoarming
            %%-------------------------------------------------------------
            BF = zeros(N_challenges,N_vox);
            for i_challenge = 1:N_challenges % challenge
                BF(i_challenge,:) = rms(estimated_sources{i_challenge},2);
            end
            %%-------------------------------------------------------------
            error = zeros(N_challenges,1);
            for i_challenge = 1:N_challenges
                %%---------------------------------------------------------
                M = BF(i_challenge,:);
                %%---------------------------------------------------------
                [~,I] = sort(M,'descend');
                order = I;
                dist = (x-circular_fascicles(i_challenge, 1)).^2 + (y-circular_fascicles(i_challenge, 2)).^2;
                dist = reshape(dist,[1 N_vox]);
                [~,I] = sort(dist,'ascend');
                cen = I(1);
                error(i_challenge) = 1e3*norm([x(cen) y(cen)]-[x(order(1,1)) y(order(1,1))]); % mm
            end
            %%-------------------------------------------------------------
            save(['SNR = ' num2str(snr) '\error_BF_' type '_' num2str(i_sec) '.mat'],'error')
            %%-------------------------------------------------------------
            
            %%-------------------------------------------------------------
            % Discriminative Field Potential
            %%-------------------------------------------------------------
            DFP = zeros(N_challenges,N_vox);
            for i_challenge = 1:N_challenges % challenge
                for j = 1:N_vox % voxels
                    DFP(i_challenge,j) = D(i_challenge,:)*hat_LFM(:,j)/(sum(hat_LFM(:,j))*sum(D(i_challenge,:)));
                    if isnan(DFP(i_challenge,j))
                        DFP(i_challenge,j) = 0;
                    end
                end
            end
            %%-------------------------------------------------------------
            error = zeros(N_challenges,1);
            for i_challenge = 1:N_challenges
                %%---------------------------------------------------------
                M = DFP(i_challenge,:);
                %%---------------------------------------------------------
                [~,I] = sort(M,'descend');
                order = I;
                dist = (x-circular_fascicles(i_challenge, 1)).^2 + (y-circular_fascicles(i_challenge, 2)).^2;
                dist = reshape(dist,[1 N_vox]);
                [~,I] = sort(dist,'ascend');
                cen = I(1);
                error(i_challenge) = 1e3*norm([x(cen) y(cen)]-[x(order(1,1)) y(order(1,1))]); % mm
            end
            %%-------------------------------------------------------------
            save(['SNR = ' num2str(snr) '\error_DFP_' type '_' num2str(i_sec) '.mat'],'error')
            %%-------------------------------------------------------------

            

            %%-------------------------------------------------------------
            % Discriminative Beamforming
            %%-------------------------------------------------------------
            DBF = zeros(N_challenges,N_vox);
            for i_challenge = 1:N_challenges % challenge
                DBF(i_challenge,:) = LFM_psd_inv*D(i_challenge,:)';
            end
            %%-------------------------------------------------------------
            error = zeros(N_challenges,1);
            for i_challenge = 1:N_challenges
                %%---------------------------------------------------------
                M = DBF(i_challenge,:);
                %%---------------------------------------------------------
                [~,I] = sort(M,'descend');
                order = I;
                dist = (x-circular_fascicles(i_challenge, 1)).^2 + (y-circular_fascicles(i_challenge, 2)).^2;
                dist = reshape(dist,[1 N_vox]);
                [~,I] = sort(dist,'ascend');
                cen = I(1);
                error(i_challenge) = 1e3*norm([x(cen) y(cen)]-[x(order(1,1)) y(order(1,1))]); % mm
            end
            %%-------------------------------------------------------------
            save(['SNR = ' num2str(snr) '\error_DBF_' type '_' num2str(i_sec) '.mat'],'error')
            %%-------------------------------------------------------------
        end
    end
end
