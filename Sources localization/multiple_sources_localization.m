clc
clear
%%-------------------------------------------------------------------------
% General info: Perform a multiple sources localization
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
% Is performed a challenge for all combination of three fascicles (at time)
% For each challenge the source signals are imposed to be a sin
% All fascicles are a sorce of brownian noise
%%-------------------------------------------------------------------------
N_sec = 1; % Define the number of section of the model to be generate
Model = 'Structural'; % 'Generic' otherwise
%%-------------------------------------------------------------------------
Data1 = 'Error 16 Hz';
Data2 = 'Error 4 Hz';
Data3 = 'Error 2 Hz';
Data4 = 'Fascicle 16 Hz';
Data5 = 'Fascicle 4 Hz';
Data6 = 'Fascicle 2 Hz';
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
    load('nerve_mod_vagus_human_generic.mat','R')
    voxels = 40; % 50 micron res
    x = linspace(-R, R, voxels);
    y = linspace(-R, R, voxels);
    z = 0;
    [x, y, ~] = meshgrid(x, y, z);
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
        l_cc = l_shaft/n;
        sel_gen = sel_gen - (x > (-1e-3*nerve_diam/2+l_cc/2) & x < (-1e-3*nerve_diam/2+l_cc/2+l_shaft) & y > (-h_shaft/2-h_as_t) & y < (-h_shaft/2-h_as_t+h_shaft));
    end
    %%---------------------------------------------------------------------

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
        %%---------------------------------------------------------------------
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
        sig = A*sin(2*pi*i_challenge*t);
        for snr = 1
            load(['SNR = ' num2str(snr) '\mod_noise_' type '_' num2str(num) '.mat'])
            errorDFP = zeros(3*N_fasc,1);
            errorBF = zeros(3*N_fasc,1);
            errorDBF = zeros(3*N_fasc,1);
            perm = nchoosek(1:N_fasc,3);
            k = 0;
            for i_fasc = 1:size(perm,1)
                permu(k+1:k+6,:) = perms(perm(i_fasc,:));
                k = k + 6;
            end
            perm = permu;
            clear permu
            N_challenges = size(perm,1);
            %%-------------------------------------------------------------
            for i_challenge = 1:N_challenges
                %%---------------------------------------------------------
                % Challenges:
                % 3 simultaneous challenges (at different frequency)
                %%---------------------------------------------------------
                N_sources = 3;
                %%---------------------------------------------------------
                % Assign which fascicle contain a source:
                first = perm(i_challenge,1);
                second = perm(i_challenge,2);
                third = perm(i_challenge,3);
                challenges = [first,second,third];
                %%---------------------------------------------------------


                %%---------------------------------------------------------
                % Rest
                %%---------------------------------------------------------
                A = 1e-6;
                rest_sig = zeros(N_vox,length(t));
                %%---------------------------------------------------------
                for i_fasc = 1:N_fasc
                    %%-----------------------------------------------------
                    % All the fascicle
                    pix = (x-circular_fascicles(i_fasc, 1)).^2 + (y-circular_fascicles(i_fasc, 2)).^2 < circular_fascicles(i_fasc, 3)^2;
                    pix = reshape(pix,[1 N_vox]);
                    pix = find(pix == 1);
                    %%-----------------------------------------------------
                    rng(s+i_fasc,'twister');
                    cn = dsp.ColoredNoise('brown','SamplesPerFrame',length(t),'NumChannels',length(pix)); % Flicker noise
                    rng(s+i_fasc,'twister');
                    rest_sig(pix,:) = mod_noise*cn()';
                end
                %%---------------------------------------------------------

                %%---------------------------------------------------------
                % Recordings
                %%---------------------------------------------------------
                f1 = 16;
                f2 = 4;
                f3 = 2;
                sig1 = A*sin(2*pi*f1.*t);
                sig2 = A*sin(2*pi*f2.*t);
                sig3 = A*sin(2*pi*f3.*t);
                source = rest_sig;
                for i_source = challenges
                    %%-----------------------------------------------------
                    % Only the central pixel
                    %%-----------------------------------------------------
                    dist = (x-circular_fascicles(i_source, 1)).^2 + (y-circular_fascicles(i_source, 2)).^2;
                    dist = reshape(dist,[1 N_vox]);
                    [~,I] = sort(dist,'ascend');
                    pix = I(1);
                    %%-----------------------------------------------------
                    if i_source == first
                        source(pix,:) = ones(length(pix),1)*sig1;
                    elseif i_source == second
                        source(pix,:) = ones(length(pix),1)*sig2;
                    else
                        source(pix,:) = ones(length(pix),1)*sig3;
                    end
                    %%-----------------------------------------------------
                end
                recording = pseudoR*source;
                %%---------------------------------------------------------
                % Source localization
                estimated_sources = LFM_psd_inv*recording;
                %%---------------------------------------------------------


                %%---------------------------------------------------------
                % Discriminative Indexes
                %%---------------------------------------------------------
                challenge = recording;
                N = length(t);
                frequency = ((1:N)-N/2)*fs/N;
                D = zeros(N_sources,N_channel);
                for i_channel = 1:N_channel
                    spectrum = abs(fftshift(fft(challenge(i_channel,:))))/N;
                    D(1,i_channel) = 2*sum(spectrum(frequency > f1-0.1 & frequency < f1+0.1));
                    D(2,i_channel) = 2*sum(spectrum(frequency > f2-0.1 & frequency < f2+0.1));
                    D(3,i_channel) = 2*sum(spectrum(frequency > f3-0.1 & frequency < f3+0.1));
                end
                %%---------------------------------------------------------
                
                %%---------------------------------------------------------
                % Beamforming
                %%---------------------------------------------------------
                BF = rms(estimated_sources,2);
                M = BF;
                %%---------------------------------------------------------
                [~,I] = sort(M,'descend');
                order = I;
                for i_source = 1:N_sources
                    dist = (x-circular_fascicles(challenges(i_source), 1)).^2 + (y-circular_fascicles(challenges(i_source), 2)).^2;
                    dist = reshape(dist,[1 N_vox]);
                    [~,I] = sort(dist,'ascend');
                    cen = I(1);
                    errorBF(N_sources*(i_challenge-1)+i_source) = error;
                end
                %%---------------------------------------------------------

                %%---------------------------------------------------------
                % Discriminative Field Potential
                %%---------------------------------------------------------
                DFP = zeros(N_sources,N_vox);
                %%---------------------------------------------------------
                for i_vox = 1:N_vox 
                    for i_source = 1:N_sources
                        DFP(i_source,i_vox) = D(i_source,:)*hat_LFM(:,i_vox)/sum(hat_LFM(:,i_vox));
                        if isnan(DFP(i_source,i_vox))
                            DFP(i_source,i_vox) = 0;
                        end
                    end
                end
                %%---------------------------------------------------------
                for i_source = 1:N_sources
                    M = DFP(i_source,:);
                    %%-----------------------------------------------------
                    [~,I] = sort(M,'descend');
                    order = I;
                    dist = (x-circular_fascicles(challenges(i_source), 1)).^2 + (y-circular_fascicles(challenges(i_source), 2)).^2;
                    dist = reshape(dist,[1 N_vox]);
                    [~,I] = sort(dist,'ascend');
                    cen = I(1);
                    errorDFP(N_sources*(i_challenge-1)+i_source) = error;
                end
                
                %%---------------------------------------------------------
                % Discriminative Beamforming
                %%---------------------------------------------------------
                 for i_source = 1:N_sources
                    DBF = LFM_psd_inv*D(i_source,:)';
                    M = DBF;
                    %%-----------------------------------------------------
                    [~,I] = sort(M,'descend');
                    order = I;
                    dist = (x-circular_fascicles(challenges(i_source), 1)).^2 + (y-circular_fascicles(challenges(i_source), 2)).^2;
                    dist = reshape(dist,[1 N_vox]);
                    [~,I] = sort(dist,'ascend');
                    cen = I(1);
                    errorDBF(N_sources*(i_challenge-1)+i_source) = error;
                 end
            end

            %%-------------------------------------------------------------
            errorBF = reshape(errorBF,[N_sources length(errorBF)/N_sources])';
            error = [{Data1} {Data2} {Data3} {Data4} {Data5} {Data6}; {errorBF(:,1)} {errorBF(:,2)} {errorBF(:,3)} {perm(:,1)} {perm(:,2)} {perm(:,3)}];
            save(['SNR = ' num2str(snr) '\error_BF_' type '_' num2str(num) '.mat'],'error')
            %%-------------------------------------------------------------
            errorDFP = reshape(errorDFP,[N_sources length(errorDFP)/N_sources])';
            error = [{Data1} {Data2} {Data3} {Data4} {Data5} {Data6}; {errorDFP(:,1)} {errorDFP(:,2)} {errorDFP(:,3)} {perm(:,1)} {perm(:,2)} {perm(:,3)}];
            save(['SNR = ' num2str(snr) '\error_DFP_' type '_' num2str(num) '.mat'],'error')
            %%-------------------------------------------------------------
            errorDBF = reshape(errorDBF,[N_sources length(errorDBF)/N_sources])';
            error = [{Data1} {Data2} {Data3} {Data4} {Data5} {Data6}; {errorDBF(:,1)} {errorDBF(:,2)} {errorDBF(:,3)} {perm(:,1)} {perm(:,2)} {perm(:,3)}];
            save(['SNR = ' num2str(snr) '\error_DBF_' type '_' num2str(num) '.mat'],'error')
            %%-------------------------------------------------------------
        end
    end
end