function [channel_estimate,channel_predict,H_test] = LTE_UL_channel_estimator(LTE_params,ChanMod,BS,rx_ref_symbols,RefSym,RefMapping,perfect_channel,sigma_n2,UE_output,nStreams,subframe_i,SNR,channel_cor,precoding_matrix,UE)
% LTE channel estimator
% [chan_output] = LTE_channel_model(BS_output, SNR)
% Author: Stefan Pratschner, spratsch@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

H_test=zeros(2400,1);
% initialize variables
N                                   = size(RefSym,1);                                       % number of scheduled subcarriers
Nslot                               = LTE_params.Nslot;                                     % number of slots per subframe
Nsub                                = LTE_params.Nsub;                                      % number of symbols per subframe
buffer_size                         = 10;                                                   % channel interpolation buffer size
lambda_f                            = LTE_params.BS_config.channel_est_frequency_smoothing; % weight of the frequency smoothing constraint
channel_estimate                    = nan(N,LTE_params.Nsub,ChanMod.nRX,nStreams);          % estimated channel
channel_predict                     = nan(N,LTE_params.Nsub,ChanMod.nRX,nStreams);          % predicted channel
channel_estimation_method           = BS.channel_estimation_method;                         % estimation method
channel_interpolation_method        = BS.channel_interpolation_method;                      % interpolation method
channel_interpolation_past_points   = BS.channel_interpolation_past_points;                 % number of previous estimates for interpolation

% reference symbol indices
switch LTE_params.CyclicPrefix
    case 'normal'   % symbols 0 to 13
        DMRS_index = 4;     % 3rd symblo per slot
        SRS_index = 14;     % 13th symbol in a subframe
    case 'extended' % smybols 0 to 11
        DMRS_index = 3;     % 2nd symbol per slot
        SRS_index = 12;     % 11th symbol per subframe
    otherwise
        error('CP not supported');
end

% prediction horizon
if LTE_params.BS_config.channel_prediction
    pred_horizon    = LTE_params.Nsub*(LTE_params.downlink_delay+1);    % prediction horizon
else
    pred_horizon    = LTE_params.Nsub;
end

% configure spatial stream separation
gamma_bar = [1,2,4,4,6,6,8,8];
if (nStreams <= 8) || (nStreams >= 0)
    LL = gamma_bar(nStreams);
else
    error('number of streams not supported');
end



% channel estimation
switch channel_estimation_method
    case 'PERFECT'
        channel_estimate = perfect_channel;     % perfect channel knowledge
        return;                                 % return perfect channel immediately, no interpolation needed
        
    case 'MMSE_2D'
        % implementation in the interpolation part below
        
    case 'MMSE'
        % generate R matrix
        H_tmp=zeros(N,Nslot,ChanMod.nRX,nStreams);
        for rr = 1:ChanMod.nRX
            % generate R matrix
            % fist slot
            R1 = zeros(N, nStreams*N);
            for lay_ind = 1:nStreams
                R1((1:N),(1+(lay_ind-1)*N:N+(lay_ind-1)*N))=diag(RefSym(:,1,lay_ind));
            end
            % second slot
            R2 = zeros(N, nStreams*N);
            for lay_ind = 1:nStreams
                R2((1:N),(1+(lay_ind-1)*N:N+(lay_ind-1)*N))=diag(RefSym(:,2,lay_ind));
            end
            
            % generate correlation matrices
            % frequency correlation
            start_sch   = find(UE_output.CH_mapping(:,1),1,'first');
            end_sch     = find(UE_output.CH_mapping(:,1),1,'last');
            
            channel_cor1 = channel_cor(start_sch:end_sch, start_sch:end_sch);
            
            % spatial correlation
            C_s = precoding_matrix'*precoding_matrix;
            
            % total correlation
            C_h = kron(C_s,channel_cor1);  
                    
            % actual MMSE estimation
            h1 = C_h * R1' * inv(R1*C_h*R1'+sigma_n2.*eye(N)) * rx_ref_symbols(:,1,rr);     % first slot
            h2 = C_h * R2' * inv(R2*C_h*R2'+sigma_n2.*eye(N)) * rx_ref_symbols(:,2,rr);     % second slot
                        
            % save results          
            H_tmp(:,1,rr,:)=reshape(h1,N,nStreams);
            H_tmp(:,2,rr,:)=reshape(h2,N,nStreams);
                    
        end % for nRx
                
    case {'LS_AV'}   % LS using orthogonal DMRS
        % individual for each slot, then average and repeat for block fading
        
        % averaging matrix
        M=kron(eye(N/LL),ones(LL,LL))./LL;
        
        % actual estimation
        for rr = 1:ChanMod.nRX
            for ll = 1:nStreams
                H_tmp(:,1,rr,ll)=(M*diag(conj(RefSym(:,1,ll)))*rx_ref_symbols(:,1,rr));
                H_tmp(:,2,rr,ll)=(M*diag(conj(RefSym(:,2,ll)))*rx_ref_symbols(:,2,rr));  %channel estimation of DMRS
            end
        end
        H_test=H_tmp(:);% 所需测试数据 2400个
       
    case 'LS_SAV'   % LS using orthogonal DMRS, sliding window, matrix implementation
        % individual for each slot, then average and repeat for block fading
        H_tmp=zeros(N,Nslot,ChanMod.nRX,nStreams);
    
        % generate matrix for that represents sum
        M_help = zeros(N);
        M_help(1:LL,1:LL) = ones(LL);
        M = zeros(N);
        for ii = 0:N-LL
            M = M + circshift(M_help, [ii ii]);
        end
        M_weight = sum(M,2);
        M = M ./ repmat(M_weight, [1 N]);   % normalize
        
        for rr = 1:ChanMod.nRX
            for ll = 1:nStreams
                H_tmp(:,1,rr,ll)=(M*diag(conj(RefSym(:,1,ll)))*rx_ref_symbols(:,1,rr));
                H_tmp(:,2,rr,ll)=(M*diag(conj(RefSym(:,2,ll)))*rx_ref_symbols(:,2,rr));
            end
        end
                
    case 'LS_QS'    % LS quadratic smoothing for one layer and one slot
        H_tmp=zeros(N,Nslot,ChanMod.nRX,nStreams);
        
        % generate D_f matrix
        D_f = -eye(N)+diag(ones(1,N-1),1);
        D_f(end,:) = [];
        
        % generate smoothing matrix
        M = pinv(eye(N)+lambda_f(nStreams).*D_f'*D_f);
        
        % actual estimation
        for rr = 1:ChanMod.nRX
            for ll = 1:nStreams
                H_tmp(:,1,rr,ll)=(M*diag(conj(RefSym(:,1,ll)))*rx_ref_symbols(:,1,rr));
                H_tmp(:,2,rr,ll)=(M*diag(conj(RefSym(:,2,ll)))*rx_ref_symbols(:,2,rr));
            end
        end         
       
    case 'LS_AQS'   % approximate quadratic inverse
        H_tmp=zeros(N,Nslot,ChanMod.nRX,nStreams);
        
        delta_n     = 10;
         
        % generate estimator matrix
        if nStreams == 1
            M_app = eye(N);
        else
            d           = 1 + 1/(2*lambda_f(nStreams));           % argument
            toep_vec    = zeros(N,1);
            for nn = 1:N
                toep_vec(nn,1) = chebyshev_2nd(d,nn);
            end
            toep_vec    = flipud(toep_vec);
            M_app       = toeplitz([toep_vec(1:delta_n);zeros(N-delta_n,1)]);
            M_app       = M_app ./ repmat(sum(M_app,2),[1,N]);  % normalization
        end
        
        
        % actual estimation
        for rr = 1:ChanMod.nRX
            for ll = 1:nStreams
                H_tmp(:,1,rr,ll)=(M_app*diag(conj(RefSym(:,1,ll)))*rx_ref_symbols(:,1,rr));
                H_tmp(:,2,rr,ll)=(M_app*diag(conj(RefSym(:,2,ll)))*rx_ref_symbols(:,2,rr));
            end
        end
            
    case{'DFT'}
        H_tmp=zeros(N,Nslot,ChanMod.nRX,nStreams);
        
        % init parameters
        D           = 1/sqrt(N) * fft(eye(N));      % DFT matrix
        win_length  = LTE_params.Ng(1);             % choose CP length as window size
        window      = [ones(win_length,1);zeros(N-2*win_length,1);ones(win_length,1)];
        
        % estimation
        for rr = 1:ChanMod.nRX
            for ll= 1:nStreams
                for mm=1:2   % slots
                    y_tmp = D'*diag(conj(RefSym(:,mm,ll)))*rx_ref_symbols(:,mm,rr);
                    H_tmp(:,mm,rr,ll) = D*(window.*y_tmp);
                end
            end
        end
        
    otherwise
        error('estimation method not supported for blockfading')
end % switch channel_estimation_method


% channel estimate interpolation
switch ChanMod.filtering
    case 'BlockFading'
        H_tmp = mean(H_tmp,2);                                  % mean over both slots to reduce noise power
        channel_estimate    = repmat(H_tmp,1,Nsub);         % repeat for blockfading
        channel_predict     = zeros(size(channel_estimate));    % no prediction in block fading
    case 'FastFading'
        switch channel_interpolation_method  % interpolation in time domain
            case {'linear', 'spline', 'pchip'}
                for rr=1:ChanMod.nRX
                    for LL = 1:nStreams
                        support = find(sum(RefMapping,1));                              % positions where we have a channel estimate
                        H_help = H_tmp(:,:,rr,LL);
                        if ~isempty(UE.ch_est_buff) && (rr <=size(UE.ch_est_buff,3)) && (LL<=size(UE.ch_est_buff,4)) && subframe_i>ceil(channel_interpolation_past_points/2)   % if we have a previous estimate...
                            for ii = 1:channel_interpolation_past_points
                                support = [support(1,1)-LTE_params.Ns, support];        % extend support to previous channel estimate
                                H_help  = [UE.ch_est_buff(:,end-ii+1,rr,LL), H_help];   % and append previous estimate
                            end
                        end
                        
                        channel_estimate_help       = interp1(support, transpose(H_help),1:pred_horizon,channel_interpolation_method, 'extrap');
                        channel_estimate(:,:,rr,LL) = transpose(channel_estimate_help(1:Nsub,:));
                        channel_predict(:,:,rr,LL)  = transpose(channel_estimate_help(pred_horizon-Nsub+1:pred_horizon,:));
                    end
                end

            case 'linLSfit' 
                for rr=1:ChanMod.nRX
                    for LL = 1:nStreams
                        support = find(sum(RefMapping,1));                              % positions where we have a channel estimate
                        H_help = H_tmp(:,:,rr,LL);
                        if ~isempty(UE.ch_est_buff) && (rr <=size(UE.ch_est_buff,3)) && (LL<=size(UE.ch_est_buff,4)) && subframe_i>ceil(channel_interpolation_past_points/2)   % if we have a previous estimate...
                            for ii = 1:channel_interpolation_past_points
                                support = [support(1,1)-LTE_params.Ns, support];        % extend support to previous channel estimate
                                H_help  = [UE.ch_est_buff(:,end-ii+1,rr,LL), H_help];   % and append previous estimate
                            end
                        end
                        
                        for kk=1:size(H_help,1)     %for each subcarrier
                            A(:,1)  = support;
                            A(:,2)  = ones(size(H_help,2),1);
                            
                            params                      = A\transpose(H_help(kk,:));
                            channel_estimate_help(:,kk) = params(1)*[1:pred_horizon] + params(2);
                            
                        end
                        channel_estimate(:,:,rr,LL) = transpose(channel_estimate_help(1:Nsub,:));
                        channel_predict             = transpose(channel_estimate_help(pred_horizon-Nsub+1:pred_horizon,:));
                    end
                end
                        
            case 'DPSS'
                for rr=1:ChanMod.nRX
                    for LL = 1:nStreams
                        support = find(sum(RefMapping,1));                              % positions where we have a channel estimate
                        H_help = H_tmp(:,:,rr,LL);
                        if ~isempty(UE.ch_est_buff) && (rr <=size(UE.ch_est_buff,3)) && (LL<=size(UE.ch_est_buff,4)) && subframe_i>ceil(channel_interpolation_past_points/2)   % if we have a previous estimate...
                            for ii = 1:channel_interpolation_past_points
                                support = [support(1,1)-LTE_params.Ns, support];        % extend support to previous channel estimate
                                H_help  = [UE.ch_est_buff(:,end-ii+1,rr,LL), H_help];   % and append previous estimate
                            end
                        end

                        % parameters
                        window_size     = size(H_help,2);               % number of DMRS within the window
                        pred_size       = 2*LTE_params.downlink_delay;  % number of DMRS to predict
                        ref_sym_dist    = LTE_params.Ns;                % distance between two DMRS
                        f_D             = LTE_params.UE_config.user_speed*LTE_params.carrier_freq_UP/LTE_params.speed_of_light;
                        
                        % generate basis
                        C = zeros((window_size + pred_size)*LTE_params.Ns);
                        for ii=0:(window_size + pred_size)*LTE_params.Ns-1
                            for nn=0:(window_size + pred_size)*LTE_params.Ns-1
                                delta = ii-nn;
                                if delta == 0
                                    C(ii+1,nn+1) = 2*f_D*mean(LTE_params.Ts);
                                else
                                    C(ii+1,nn+1) = sin(2*pi*delta*f_D*mean(LTE_params.Ts))/(pi*delta);
                                end
                            end
                        end
                        
                        [V,D_eig] = eig(C);         % basis are eigenvectors
                        V = fliplr(V);              % sort eigenvectors
                        D_eig = rot90(D_eig,2);     % sort eigenvalues
                        F = V(:,1:window_size).';   % basis
                        
                        % correlation matrix
                        G = zeros(window_size);
                        for uuu = (support+(window_size-2)*LTE_params.Ns)
                            G = G + F(:,uuu)*F(:,uuu)';
                        end
                        G_inv = pinv(G);
                        
                        channel_estimate_help = zeros(N,(window_size + pred_size)*LTE_params.Ns);
                        for kk=1:N
                            gamma_hat = zeros(window_size,1);
                            for ooo = 1:window_size
                                gamma_hat = gamma_hat + H_help(kk,ooo)*conj(F(:,support(ooo)+(window_size-2)*LTE_params.Ns));
                            end
                            gamma_hat = G_inv * gamma_hat;

                            % interpolation
                            for tt =1:window_size
                                channel_estimate_help(kk,:) = channel_estimate_help(kk,:) + transpose(V(:,tt) * gamma_hat(tt));
                            end
                        end
                        
                        channel_estimate(:,:,rr,LL) = channel_estimate_help(:,(size(H_help,2)-2)*LTE_params.Ns+1:size(H_help,2)*LTE_params.Ns);
                        channel_predict(:,:,rr,LL)  = channel_estimate_help(:,end-14+1:end);
                    end
                end          
                
            case 'MMSE_2D'
                % generate reference symbol matrix
                RefSym_help = RefSym;
                if ~isempty(UE.refsym_buff) && (subframe_i>ceil(channel_interpolation_past_points/2))   % if there are previous received symbols
                    for ii = 1:channel_interpolation_past_points
                        RefSym_help = [UE.refsym_buff(:,end-ii+1,:) ,RefSym_help];
                    end
                end

                R = zeros(N,nStreams*N,size(RefSym_help,2));
                for ss =1:size(RefSym_help,2)
                    for lay_ind = 1:nStreams
                        R((1:N),(1+(lay_ind-1)*N:N+(lay_ind-1)*N),ss)=diag(RefSym_help(:,ss,lay_ind));
                    end
                end
                
                R_tot = blkdiag(R(:,:,1),R(:,:,2));
                for ss = 1:size(RefSym_help,2)-2
                    R_tot = blkdiag(R_tot,R(:,:,ss+2));
                end
                
                for rr=1:ChanMod.nRX
                    % consider previous received symbols
                    support = find(sum(RefMapping,1));
                    y_help = reshape(rx_ref_symbols(:,:,rr),[],1);
                    if ~isempty(UE.rx_refsym_buff) && (subframe_i>ceil(channel_interpolation_past_points/2))   % if there are previous received symbols
                        for ii = 1:channel_interpolation_past_points
                            support = [support(1,1)-LTE_params.Ns, support];        % extend support to previous channel estimate
                            y_help = [UE.rx_refsym_buff(:,end-ii+1,rr) ;y_help];    % prepend previos received symbols
                        end
                    end
                                        
                    % parameters
                    window_size     = size(y_help,1)/N;             % number of DMRS within the window
                    pred_size       = 2*LTE_params.downlink_delay;  % number of DMRS to predict
                    ref_sym_dist    = LTE_params.Ns;                % distance between two DMRS
                    f_D             = LTE_params.UE_config.user_speed*LTE_params.carrier_freq_UP/LTE_params.speed_of_light;

                    % frequency correlation
                    start_sch   = find(UE_output.CH_mapping(:,1),1,'first');
                    end_sch     = find(UE_output.CH_mapping(:,1),1,'last');

                    channel_cor1 = channel_cor(start_sch:end_sch, start_sch:end_sch);
                    
                    % no spatial correlation
                    C_h = kron(eye(nStreams),channel_cor1);

                    % time correlation
                    % generate autocorrelation matrix
                    R_HH = zeros(window_size);
                    for uu = 1:window_size
                        for vv = 1:window_size
                            R_HH(uu,vv) = besselj(0,2*pi*abs(uu-vv)*ref_sym_dist*f_D*mean(LTE_params.Ts));
                        end
                    end 

                    % generate crosscorrelation vector
                    r_Hn = zeros(window_size,(window_size+pred_size)*ref_sym_dist);
                    for tt = 1:(window_size+pred_size)*ref_sym_dist
                        for ii = 1:window_size
                            r_Hn(ii,tt) = besselj(0,2*pi*abs((ii-1)*ref_sym_dist-(tt-4))*f_D*mean(LTE_params.Ts));
                        end
                    end

                    % calculate vectorized correlation matrices
                    C_auto_vec  = R_tot*kron(R_HH,C_h)*R_tot' + sigma_n2*eye(window_size*N);
                    C_cross_vec = kron(r_Hn.',C_h)*R_tot';          

                    % construct receive vector
%                     y_stack = reshape(rx_ref_symbols(:,:,rr),[],1);
                    
                    % 2D MMSE interpolation
                    channel_estimate_vec    = C_cross_vec * inv(C_auto_vec) * y_help;

                    channel_estimate_help(:,:,rr,:) = permute(reshape(channel_estimate_vec, N, nStreams, []),[1,3,2]);
                    channel_estimate(:,:,rr,:)      = channel_estimate_help(:,(window_size-2)*LTE_params.Ns+1:window_size*LTE_params.Ns,rr,:);
                    channel_predict(:,:,rr,:)       = channel_estimate_help(:,end-14+1:end,rr,:);
                end
                
                % cunstruct H_tmp
                H_tmp = channel_estimate(:,[DMRS_index,DMRS_index+LTE_params.Ns],:,:);
                
            case 'flat' % repeat channel estimate for one slot
                for rr=1:ChanMod.nRX
                    for LL = 1:nStreams
                        H_help = H_tmp(:,:,rr,LL);
                        channel_estimate(:,:,rr,LL) = kron(H_help, ones(1,LTE_params.Ns));
                        channel_predict(:,:,rr,LL)  = repmat(H_help(:,2), 1, LTE_params.Nsub);
                        
                    end
                end
                
            case 'SVM'
                for rr=1:ChanMod.nRX
                    for LL = 1:nStreams
                        support = find(sum(RefMapping,1));                          % positions where we have a channel estimate
                        H_help = H_tmp(:,:,rr,LL);
                       if ~isempty(UE.ch_est_buff) && (rr <=size(UE.ch_est_buff,3)) && (LL<=size(UE.ch_est_buff,4)) && subframe_i>ceil(channel_interpolation_past_points/2)   % if we have a previous estimate...
                            for ii = 1:channel_interpolation_past_points
                                support = [support(1,1)-LTE_params.Ns, support];        % extend support to previous channel estimate
                                H_help  = [UE.ch_est_buff(:,end-ii+1,rr,LL), H_help];   % and append previous estimate
                            end
                        end
                        
                        [channel_estimate,channel_predict]=svr_pred(H_help,support,SNR,subframe_i,rr,LL,pred_horizon,Nsub);
                    end
                end
            otherwise
                error('interpolation method not supported')
        end
        
        if channel_interpolation_past_points    % in case previous points are considered
            % clear all buffers when channel dimension (nLayers) has changed
            if(size(UE.ch_est_buff,4)~=size(H_tmp,4))                                           % if the channel dimension (nLayers) changed
                UE.ch_est_buff      = [];                                                       % clear the buffer
                UE.rx_refsym_buff   = [];                                                       % clear the buffer
                UE.refsym_buff      = [];                                                       % clear the buffer
            end
            
            % save channel estimate to genie for interpolation of the next subframe
            if(size(UE.ch_est_buff,2)<buffer_size)                                              % as long as the buffer is not full
                UE.ch_est_buff = [UE.ch_est_buff, H_tmp];                                       % append the current estimate
            else                                                                                % if the buffer is full
                UE.ch_est_buff = [UE.ch_est_buff(:,3:buffer_size,:,:),H_tmp];                   % append the current estimate, and drop the oldest estimate
            end

            % save received reference symbols
            if(size(UE.rx_refsym_buff,2)<buffer_size)                                           % as long as the buffer is not full
                UE.rx_refsym_buff = [UE.rx_refsym_buff, rx_ref_symbols ];                       % append the current received symbols
            else                                                                                % if the buffer is full
                UE.rx_refsym_buff = [UE.rx_refsym_buff(:,3:buffer_size,:),rx_ref_symbols];      % append the received symbols, and drop the oldest
            end

            % save reference symbols
            if(size(UE.refsym_buff,2)<buffer_size)                                              % as long as the buffer is not full
                UE.refsym_buff = [UE.refsym_buff, RefSym ];                                     % append the current reference symbols
            else                                                                                % if the buffer is full
                UE.refsym_buff = [UE.refsym_buff(:,3:buffer_size,:),RefSym];                    % append the received symbols, and drop the oldest
            end
        end

    otherwise
        error('selected filtering type not supported');
end

% save channel estimate to genie for interpolation of the next subframe
% UE_genie.channel_estimate = H_tmp;
