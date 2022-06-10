clc
clear
%%-------------------------------------------------------------------------
% General info: Perform a multiple biophysical sources localization
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
% Analysis are done with a dt = 0.025 ms equally to fs = 40 kHz
%%-------------------------------------------------------------------------
% Challenges:
% For each afferent fascicle (at time) signal is imposed to be a spike
% train modeled as an inhomogeneous Poisson process, with an imposed 
% firing rate between 2.5 and 25 Hz, modulated by physiological  
% parameters (BP and RESP)
% Others fascicles activity is modeled as a homogeneous Poisson process
% with an imposed firing rate between 2.5 and 25 Hz
%%-------------------------------------------------------------------------
N_sec = 1; % Define the number of section of the model to be generate
Model = 'Structural'; % 'Generic' otherwise
%%-------------------------------------------------------------------------
Data1 = 'Error BP';
Data2 = 'Error RESP';
Data3 = 'SNR BP';
Data4 = 'SNR RESP';
Data5 = 'Fascicle BP';
Data6 = 'Fascicle RESP';
Legend = {Data1 Data2 Data3 Data4 Data5 Data6};
%%-------------------------------------------------------------------------
for e = 1:2
    switch e
        case 1
            type = 'TIME';
        case 2
            type = 'CUFF';
    end
    load('nerve_mod_vagus_human_generic.mat','R')
    voxels = 40; % 50 micron res
    x = linspace(-R, R, voxels);
    y = linspace(-R, R, voxels);
    z = 0;
    [x, y, z] = meshgrid(x, y, z);
    xq = x(:);
    yq = y(:);
    zq = z(:);
    sel = x.^2 + y.^2 < R^2;
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
        sel = sel - (x > (-1e-3*nerve_diam/2+l_cc/2) & x < (-1e-3*nerve_diam/2+l_cc/2+l_shaft) & y > (-h_shaft/2-h_as_t) & y < (-h_shaft/2-h_as_t+h_shaft));
    end
    %%---------------------------------------------------------------------

    %%---------------------------------------------------------------------
    load('syn_bp.mat')    
    load('syn_resp.mat')
    %%---------------------------------------------------------------------
    % [crossCorr,lags] = xcorr(bp,resp,3/dt,'unbiased');
    % [~,I] = max(crossCorr);
    % tau = lags(I)+(60/75)/dt;
    tau = 107794;
    %%---------------------------------------------------------------------
    % Allign BP and RESP to the Respiratory Sinus Arrhythmia
    %%---------------------------------------------------------------------
    bp = bp(1+10*(1/dt)-tau:end);
    bp = bp(1:end-tau);
    %%---------------------------------------------------------------------
    resp = resp(1+10*(1/dt)-tau:end);
    resp = resp(tau+1:end);  
    %%---------------------------------------------------------------------
    for i_sec = 1:N_sec % current number of nerve model
        %%-----------------------------------------------------------------
        % Integrate structural information
        %%-----------------------------------------------------------------
        load(['nerve_mod_vagus_human_' num2str(i_sec) '.mat'],'circular_fascicles')
        N_fasc = size(circular_fascicles,1);
        if strcmp(type,'TIME')
            circular_fascicles = circular_fascicles_TIME;
        end
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
        %%---------------------------------------------------------------------
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
        load(['nerve_mod_vagus_human_' num2str(i_sec) '_' type '_rec'],'motor_fasc')
        %%-----------------------------------------------------------------
        N_ch = size(LFM_psd_inv,2);
        N_vox = size(LFM_psd_inv,1);
        %%-----------------------------------------------------------------
        % Signal processing parameters 
        %%-----------------------------------------------------------------
        % Enhanced BP
        fcb = 10;
        mawb = 0.15;
        hpb = 1000;
        lpb = 6000;
        %%-----------------------------------------------------------------
        % Enhanced RESP
        fcr = 4;
        mawr = 0.15;
        hpr = 1000;
        lpr = 6000;
        %%-----------------------------------------------------------------

        %%-----------------------------------------------------------------
        fasc = 1:N_fasc;
        fasc = setdiff(fasc,motor_fasc);
        perm = nchoosek(fasc,2);
        k = 0;
        for i_fasc = 1:size(perm,1)
            permu(k+1:k+2,:) = perms(perm(i_fasc,:));
            k = k + 2;
        end
        perm = permu;
        clear permu
        N_challenges = size(perm,1);
        %%-----------------------------------------------------------------
        SNR = zeros(N_challenges,2);
        D = zeros(N_fasc,N_ch);
        weigth = zeros(N_fasc,N_ch);
        err = zeros(N_fasc,N_ch);
        BF = zeros(N_challenges,N_vox);
        error_bf_bp = zeros(N_challenges,1);
        error_bf_resp = zeros(N_challenges,1);
        DFP = zeros(N_challenges,N_vox);
        error_dfp_bp = zeros(N_challenges,1);
        error_dfp_resp = zeros(N_challenges,1);
        DBF = zeros(N_challenges,N_vox);
        error_dbf_bp = zeros(N_challenges,1);
        error_dbf_resp = zeros(N_challenges,1);
        %%-----------------------------------------------------------------
        for i_challenge = 1:N_challenges
            load(['recordings_' type '_human_vagus_' num2str(i_sec) '_' num2str(i_challenge) '.mat'],'data','bp_fasc','resp_fasc','noise_fasc_sig','bp_fasc_sig','resp_fasc_sig')
            sig = data;
            SNRBP = mean(rms(bp_fasc_sig,2))/mean(rms(noise_fasc_sig,2)+rms(resp_fasc_sig,2));
            SNRRESP = mean(rms(resp_fasc_sig,2))/mean(rms(noise_fasc_sig,2)+rms(bp_fasc_sig,2));
            SNR(i_challenge,:) = [SNRBP SNRRESP];        
            %%-------------------------------------------------------------
            % Source localization
            estimated_sources = LFM_psd_inv*sig(:,1:(1/dt));
            %%-------------------------------------------------------------
            

            %%-------------------------------------------------------------
            % Pre-processing the neural signals for the neural power
            % envelope (of the blood pressure)
            %%-------------------------------------------------------------
            [b,a] = ellip(4,0.1,50,2*[hpb lpb]/fs,'bandpass');
            data = filtfilt(b,a,data')';
            rms_sig = sqrt(movmean(data.*data,mawb*(1/dt),2));
            [blow,alow] = butter(4,fcb/(fs/2));
            rms_sig = filtfilt(blow,alow,rms_sig')';
            %%-------------------------------------------------------------


            %%-------------------------------------------------------------
            % Discriminative Indexes
            %%-------------------------------------------------------------
            % Linear model: Xeng = D*H + e with H = [1 Xbp Xresp]
            %%-------------------------------------------------------------
            H = [ones(trec*(1/dt),1) bp resp];
            % removing filter artifacts
            H = H(1/dt+1:end-1/dt,:);
            rms_sig = rms_sig(:,1/dt+1:end-1/dt);
            rms_sig = rms_sig';
            beta = pinv(H)*rms_sig;
            D(i_challenge,:) = beta(2,:);
            %%-------------------------------------------------------------

            %%-------------------------------------------------------------
            % Beamforming
            %%-------------------------------------------------------------
            BF(i_challenge,:) = rms(estimated_sources,2);
            %%-------------------------------------------------------------
            M = BF(i_challenge,:);
            [~,I] = sort(M,'descend');
            order = I;
            dist = (x-circular_fascicles(bp_fasc, 1)).^2 + (y-circular_fascicles(bp_fasc, 2)).^2;
            dist = reshape(dist,[1 N_vox]);
            [~,I] = sort(dist,'ascend');
            cen = I(1);
            error_bf_bp(i_challenge) = 1e3*norm([x(cen) y(cen)]-[x(order(1,1)) y(order(1,1))]); % mm
            %%-------------------------------------------------------------
            
            %%-------------------------------------------------------------
            % Discriminative Field Potential
            %%-------------------------------------------------------------
            for i_fib = 1:N_vox % voxels
                DFP(i_challenge,i_fib) = D(i_challenge,:)*hat_LFM(:,i_fib)/sum(hat_LFM(:,i_fib));
                if isnan(DFP(i_challenge,i_fib))
                    DFP(i_challenge,i_fib) = 0;
                end
            end
            %%-------------------------------------------------------------
            M = DFP(i_challenge,:);
            %%-------------------------------------------------------------
            [~,I] = sort(M,'descend');
            order = I;
            dist = (x-circular_fascicles(bp_fasc, 1)).^2 + (y-circular_fascicles(bp_fasc, 2)).^2;
            dist = reshape(dist,[1 N_vox]);
            [~,I] = sort(dist,'ascend');
            cen = I(1);
            error_dfp_bp(i_challenge) = 1e3*norm([x(cen) y(cen)]-[x(order(1,1)) y(order(1,1))]); % mm
            %%-------------------------------------------------------------
            
            %%-------------------------------------------------------------
            % Discriminative Beamforming
            %%-------------------------------------------------------------
            DBF(i_challenge,:) = LFM_psd_inv*D(i_challenge,:)';
            %%-------------------------------------------------------------
            M = DBF(i_challenge,:);
            [~,I] = sort(M,'descend');
            order = I;
            dist = (x-circular_fascicles(bp_fasc, 1)).^2 + (y-circular_fascicles(bp_fasc, 2)).^2;
            dist = reshape(dist,[1 N_vox]);
            [~,I] = sort(dist,'ascend');
            cen = I(1);
            error_dbf_bp(i_challenge) = 1e3*norm([x(cen) y(cen)]-[x(order(1,1)) y(order(1,1))]); % mm
            %%-------------------------------------------------------------

            %%-------------------------------------------------------------
            % Pre-processing the neural signals for the neural power
            % envelope (of the respiration)
            %%-------------------------------------------------------------
            [b,a] = ellip(4,0.1,50,2*[hpr lpr]/fs,'bandpass');
            data = filtfilt(b,a,sig')';
            rms_sig = sqrt(movmean(data.*data,mawr*(1/dt),2));
            [blow,alow] = butter(4,fcr/(fs/2));
            rms_sig = filtfilt(blow,alow,rms_sig')';
            %%-------------------------------------------------------------


            %%-------------------------------------------------------------
            % Discriminative Indexes
            %%-------------------------------------------------------------
            % Linear model: Xeng = D*H + e with H = [1 Xbp Xresp]
            %%-------------------------------------------------------------
            H = [ones(trec*(1/dt),1) bp resp];
            % removing filter artifacts
            H = H(1/dt+1:end-1/dt,:);
            rms_sig = rms_sig(:,1/dt+1:end-1/dt);
            rms_sig = rms_sig';
            beta = pinv(H)*rms_sig;
            D(i_challenge,:) = beta(3,:);
            %%-------------------------------------------------------------

            %%-------------------------------------------------------------
            % Beamforming
            %%-------------------------------------------------------------
            BF(i_challenge,:) = rms(estimated_sources,2);
            %%-------------------------------------------------------------
            M = BF(i_challenge,:);
            [~,I] = sort(M,'descend');
            order = I;
            dist = (x-circular_fascicles(resp_fasc, 1)).^2 + (y-circular_fascicles(resp_fasc, 2)).^2;
            dist = reshape(dist,[1 N_vox]);
            [~,I] = sort(dist,'ascend');
            cen = I(1);
            error_bf_resp(i_challenge) = 1e3*norm([x(cen) y(cen)]-[x(order(1,1)) y(order(1,1))]); % mm
            %%-------------------------------------------------------------

            %%-------------------------------------------------------------
            % Discriminative Field Potential
            %%-------------------------------------------------------------
            for i_fib = 1:N_vox % voxels
                DFP(i_challenge,i_fib) = D(i_challenge,:)*hat_LFM(:,i_fib)/sum(hat_LFM(:,i_fib));
                if isnan(DFP(i_challenge,i_fib))
                    DFP(i_challenge,i_fib) = 0;
                end
            end
            %%-------------------------------------------------------------
            M = DFP(i_challenge,:);
            %%-------------------------------------------------------------
            [~,I] = sort(M,'descend');
            order = I;
            dist = (x-circular_fascicles(resp_fasc, 1)).^2 + (y-circular_fascicles(resp_fasc, 2)).^2;
            dist = reshape(dist,[1 N_vox]);
            [~,I] = sort(dist,'ascend');
            cen = I(1);
            error_dfp_resp(i_challenge) = 1e3*norm([x(cen) y(cen)]-[x(order(1,1)) y(order(1,1))]); % mm
            %%-------------------------------------------------------------

            

            %%-------------------------------------------------------------
            % Discriminative Beamforming
            %%-------------------------------------------------------------
            DBF(i_challenge,:) = LFM_psd_inv*D(i_challenge,:)';
            %%-------------------------------------------------------------
            M = DBF(i_challenge,:);
            [~,I] = sort(M,'descend');
            order = I;
            dist = (x-circular_fascicles(resp_fasc, 1)).^2 + (y-circular_fascicles(resp_fasc, 2)).^2;
            dist = reshape(dist,[1 N_vox]);
            [~,I] = sort(dist,'ascend');
            cen = I(1);
            error_dbf_resp(i_challenge) = 1e3*norm([x(cen) y(cen)]-[x(order(1,1)) y(order(1,1))]); % mm
            %%-------------------------------------------------------------

        end
        %%-----------------------------------------------------------------
        error = [{Data1} {Data2} {Data3} {Data4} {Data5} {Data6}; {error_bf_bp} {error_bf_resp} {SNR(:,1)} {SNR(:,2)} {perm(:,1)} {perm(:,2)}];
        save(['error_BF_' type '_' num2str(i_sec) '.mat'],'error')
        %%-----------------------------------------------------------------
        error = [{Data1} {Data2} {Data3} {Data4} {Data5} {Data6}; {error_dfp_bp} {error_dfp_resp} {SNR(:,1)} {SNR(:,2)} {perm(:,1)} {perm(:,2)}];
        save(['error_DFP_' type '_' num2str(i_sec) '.mat'],'error')
        %%---------------------------------------------------------------------
        error = [{Data1} {Data2} {Data3} {Data4} {Data5} {Data6}; {error_dbf_bp} {error_dbf_resp} {SNR(:,1)} {SNR(:,2)} {perm(:,1)} {perm(:,2)}];
        save(['error_DBF_' type '_' num2str(i_sec) '.mat'],'error')
        %%-----------------------------------------------------------------
    end
end
