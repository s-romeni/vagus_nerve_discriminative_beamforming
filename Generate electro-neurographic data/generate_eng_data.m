clc
clear
%%-------------------------------------------------------------------------
% General info: generate synthetic ENG signals from the left human vagus
% nerve using a ihnomogenous Poisson process
%%-------------------------------------------------------------------------
% Along the scripts several parameters have been setted based on the
% following bibliographic source:
% M. M. Ottaviani, L. Wright, T. Dawood, and V. G. Macefield, “In vivo
% recordings from the human vagus nerve using ultrasound‐guided
% microneurography,” J Physiol, vol. 598, no. 17, pp. 3569–3576, Sep. 2020,
% doi: 10.1113/JP280077.
%%-------------------------------------------------------------------------
% Authors: 
%%-------------------------------------------------------------------------
% Andrea Pitzus @TNE, SSSA // @MeDSP, UniCa & Simone Romeni @TNE, EPFL
%%-------------------------------------------------------------------------
% Analysis are done with a dt = 0.025 ms equally to fs = 40 kHz
%%-------------------------------------------------------------------------
N_sec = 1; % Define the number of section of the model to be generate
%%-------------------------------------------------------------------------

%%-------------------------------------------------------------------------
for e = 1:2
    switch e
        case 1
            type = 'TIME';
        case 2
            type = 'CUFF';
    end
    %%---------------------------------------------------------------------
    % Challenges:
    % For each afferent fascicle (at time) signal is imposed to be a spike
    % train modeled as an inhomogeneous Poisson process, with an imposed 
    % firing rate between 2.5 and 25 Hz, modulated by physiological  
    % parameters (BP and RESP)
    % Others fascicles activity is modeled as a homogeneous Poisson process
    % with an imposed firing rate between 2.5 and 25 Hz
    %%---------------------------------------------------------------------
    load('syn_bp.mat')
    load('syn_resp.mat')
    %%---------------------------------------------------------------------
    dt = t(2)-t(1);
    fs = 1/dt;
    trec = 30;
    t = t(1:trec*(1/dt));
    %%---------------------------------------------------------------------
    % [crossCorr,lags] = xcorr(bp,resp,3/dt,'unbiased');
    % [~,I] = max(crossCorr);
    % tau = lags(I)+(60/75)/dt;
    tau = 107794;
    %%---------------------------------------------------------------------
    % Prepare the time-dependent rate function
    %%---------------------------------------------------------------------
    bp = bp(1+10*(1/dt)-tau:end);
    bp = bp(1:end-tau);
    bp = bp';
    lambdab = (bp-min(bp))/(max(bp)-min(bp));
    lambdab = (22.5)*lambdab+2.5;
    %%---------------------------------------------------------------------
    resp = resp(1+10*(1/dt)-tau:end);
    resp = resp(tau+1:end);
    resp = resp';
    lambdar = (resp-min(resp))/(max(resp)-min(resp));
    lambdar = (22.5)*lambdar+2.5;
    %%---------------------------------------------------------------------
    lam = 25;
    %%---------------------------------------------------------------------

    %%---------------------------------------------------------------------
    for i_sec = 1:N_sec % current number of nerve model
        %%-----------------------------------------------------------------
        load(['nerve_mod_vagus_human_' num2str(i_sec) '.mat'],'circular_fascicles','circular_fascicles_TIME')
        if strcmp(type,'TIME')
            circular_fascicles = circular_fascicles_TIME;
        end
        N_fasc = size(circular_fascicles,1);
        %%-----------------------------------------------------------------
        load(['nerve_mod_vagus_human_' num2str(i_sec) '_' type '_rec'],'nerve','motor_fasc')
        %%-----------------------------------------------------------------

        
        %%-----------------------------------------------------------------
        fasc = 1:N_fasc;
        fasc = setdiff(fasc,motor_fasc);
        perm = nchoosek(fasc,2);
        k = 0;
        for i = 1:size(perm,1)
            permu(k+1:k+2,:) = perms(perm(i,:));
            k = k + 2;
        end
        perm = permu;
        clear permu
        %%-----------------------------------------------------------------

        %%-----------------------------------------------------------------
        N_challenge = size(perm,1);
        %%-----------------------------------------------------------------
        for i_chal = 1:N_challenge
            sig = zeros(14,length(t));
            bp_fasc_sig = zeros(14,length(t));
            resp_fasc_sig = zeros(14,length(t));
            noise_fasc_sig = zeros(14,length(t));
            bp_fasc = perm(i_chal,1);
            resp_fasc = perm(i_chal,2);
            if ~exist(['recordings_' type '_human_vagus_' num2str(i_sec) '_' num2str(i_chal) '.mat'], 'file')
                for i_fasc = 1:N_fasc
                    N_fibers = length(nerve.fascicles{i_fasc});
                    if i_fasc == bp_fasc % blood pressure modulation
                        for i_fib = 1:N_fibers
                            u = rand(1,ceil(1.5*trec*lam));
                            events = cumsum(-(1/lam)*log(u));
                            events = events(events<trec);
                            N = length(events);
                            events = events(rand(1,N)<lambdab(ceil(events/dt))/lam);
                            events = sort(events);
                            delay = max(abs(nerve.fascicles{i_fasc}(i_fib).data(1,:)));
                            %%---------------------------------------------
                            % Check refractory period
                            isi = events(2:end)-events(1:end-1);
                            rem = find(isi < 1.5e-3);
                            events(rem+1) = [];
                            events(events > t(end-1e2)) = [];
                            events = events-delay;
                            %%---------------------------------------------
                            temp = zeros(14,length(t));
                            for l = 1:length(events)
                                time_wind = floor(events(l)/dt)+1:floor(events(l)/dt)+size(nerve.fascicles{i_fasc}(i_fib).data,2);
                                temp(:,time_wind) = temp(:,time_wind) + nerve.fascicles{i_fasc}(i_fib).data;
                            end
                            sig = sig + temp;
                            bp_fasc_sig = bp_fasc_sig + temp;
                        end
                    elseif i_fasc == resp_fasc % respiration modulation
                        for i_fib = 1:N_fibers
                            u = rand(1,ceil(1.5*trec*lam));
                            events = cumsum(-(1/lam)*log(u));
                            events = events(events<trec);
                            N = length(events);
                            events = events(rand(1,N)<lambdar(ceil(events/dt))/lam);
                            events = sort(events);
                            delay = max(abs(nerve.fascicles{i_fasc}(i_fib).data(1,:)));
                            %%---------------------------------------------
                            % Check refractory period
                            isi = events(2:end)-events(1:end-1);
                            rem = find(isi < 1.5e-3);
                            events(rem+1) = [];
                            events(events > t(end-1e2)) = [];
                            events = events-delay;
                            %%---------------------------------------------
                            temp = zeros(14,length(t));
                            for l = 1:length(events)
                                time_wind = floor(events(l)/dt)+1:floor(events(l)/dt)+size(nerve.fascicles{i_fasc}(i_fib).data,2);
                                temp(:,time_wind) = temp(:,time_wind) + nerve.fascicles{i_fasc}(i_fib).data;
                            end
                            sig = sig + temp;
                            resp_fasc_sig = resp_fasc_sig + temp;
                        end
                    else
                        for i_fib = 1:N_fibers
                            lambda = (25-2.5)*rand(1)+2.5;
                            r = poissrnd(trec*lambda);
                            events = trec*abs(rand(floor(r),1));
                            events = sort(events);
                            delay = max(abs(nerve.fascicles{i_fasc}(i_fib).data(1,:)));
                            %%---------------------------------------------
                            % Check refractory period
                            isi = events(2:end)-events(1:end-1);
                            rem = find(isi < 1.5e-3);
                            events(rem+1) = [];
                            events(events > t(end-1e2)) = [];
                            events = events-delay;
                            %%---------------------------------------------
                            temp = zeros(14,length(t));
                            for l = 1:length(events)
                                time_wind = floor(events(l)/dt)+1:floor(events(l)/dt)+size(nerve.fascicles{i_fasc}(i_fib).data,2);
                                temp(:,time_wind) = temp(:,time_wind) + nerve.fascicles{i_fasc}(i_fib).data;
                            end
                            sig = sig+temp;
                            noise_fasc_sig = noise_fasc_sig + temp;
                        end
                    end
                end
                % SNR_BP = mean(rms(bp_fasc_sig,2))/mean(rms(noise_fasc_sig,2)+rms(resp_fasc_sig,2));
                % SNR_RESP = mean(rms(resp_fasc_sig,2))/mean(rms(noise_fasc_sig,2)+rms(bp_fasc_sig,2));
                % SNR(i_chal,:) = [SNR_BP SNR_RESP];
                %%---------------------------------------------------------
                data = sig;
                save(['recordings_' type '_human_vagus_' num2str(i_sec) '_' num2str(i_chal) '.mat'],'data','bp_fasc','resp_fasc','noise_fasc_sig','bp_fasc_sig','resp_fasc_sig')
            else
                %%---------------------------------------------------------
                load(['recordings_' type '_human_vagus_' num2str(i_sec) '_' num2str(i_chal) '.mat'],'data','bp_fasc','resp_fasc','noise_fasc_sig','bp_fasc_sig','resp_fasc_sig')
                sig = data;
                % SNR_BP = mean(rms(bp_fasc_sig,2))/mean(rms(noise_fasc_sig,2)+rms(resp_fasc_sig,2));
                % SNR_RESP = mean(rms(resp_fasc_sig,2))/mean(rms(noise_fasc_sig,2)+rms(bp_fasc_sig,2));
                % SNR(i_chal,:) = [SNR_BP SNR_RESP];
                %%---------------------------------------------------------
            end
        end
    end
end
