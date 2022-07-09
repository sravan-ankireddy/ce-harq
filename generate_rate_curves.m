% Fixing the range for reproducibility 
rng("default");

% Simulation settings : Defaults
if ~exist('global_settings','var')
    % Code parameters
    nPRB = 20; % Vary this to change the code length
    max_iter = 6; % default is 8 in MATLAB
    num_blocks = 10e3;
end

% Simulation params
SNRdB_step = 0.2;
SNRdB_low = -10;
SNRdB_high = 20;
SNRdB_vec = SNRdB_low:SNRdB_step:SNRdB_high;
num_snr = length(SNRdB_vec);

tcr_list = 0.1:0.1:0.9;

% data array : Store SNR,Rate,BER
data_array = zeros(length(tcr_list),num_snr);

parpool_size = min(64,feature('numcores'));
if (parpool_size ~= 8)
    parpool(parpool_size);
end
tic;
for it = 1:length(tcr_list)
    targetCodeRate = tcr_list(it);
    
    % use nrTBS to get K,N
    modulation = 'pi/2-BPSK';
    nlayers = 1;
    
    NREPerPRB = 12*4; % For URLLC, 2-7 is the typical choice
    tbs = nrTBS(modulation,nlayers,nPRB,NREPerPRB,targetCodeRate);
    
    % Use generated code params
    K = tbs;
    N = nPRB*NREPerPRB;
    R = K/N;
    
    % base graph selection based on rate and K
    bgn = bgn_select(K,R);
    
    rv = 0;

    M = 2;
    
    BER_vec = zeros(size(SNRdB_vec));
    BLER_vec = zeros(size(SNRdB_vec));

    parfor i_s = 1:length(SNRdB_vec)
        SNRdB = SNRdB_vec(i_s);
        noiseVar = 1./(10.^(SNRdB/10));
    
        BER = 0; BLER = 0;
    
        % Print status
        fprintf('Status %0.2f %% done \n', round(i_s/num_snr*100));
        
        for i_n = 1:num_blocks
            % Generate random message
            data = randi([0 1], K, 1);
        
            % Encoding and Rate matching
            bgn = bgn_select(K,R);
            dataIn = nrldpc_enc(data, R, modulation, rv, bgn, nlayers);
    
            % Reshape data into binary k-tuples, k = log2(M), and convert to
            % integers
            k = log2(M);
            txData = reshape(dataIn,length(dataIn)/k,k);
            txDataSym = bi2de(txData);
    
            % QAM Modulation
            txSig = qammod(txDataSym,M);
                
            % Pass through AWGN channel
            rxSig = awgn(txSig,SNRdB);
    
            % QAM Demod
            rxLLR = qamdemod(rxSig,M,'OutputType','LLR');
                    
            % Decoding and Rate recovery
            bgn = bgn_select(K,R);
            [data_est, crc_chk] = nrldpc_dec(rxLLR, R, modulation, K, max_iter, rv, nlayers, bgn);
    
            % Check for error stats
            num_err = sum(mod(data+double(data_est),2));
    
            BER = BER + num_err;
    
            if (num_err > 0)
                BLER = BLER + 1;
            end
        end
    
        BER_vec(i_s) = BER/(K*num_blocks);
        BLER_vec(i_s) = BLER/(num_blocks);

    end

    % Store the data
    data_array(it,:) = BER_vec;
    
    data = data_array;
    data = [transpose(tcr_list) data];
    data = [[0 SNRdB_vec]; data];
    save('results_rate_curves/BER_LDPC_960.mat','data');
    xlabel('SNR');
    figure(1);
    semilogy(SNRdB_vec,BER_vec,'-o');
    hold on;
    f = semilogy(SNRdB_vec,BLER_vec,'-d');
    
end

toc;
if (parpool_size ~= 8)
    delete(gcp('nocreate'));
end

xlabel('SNR');
ylabel('Error Rates');
legend('BER R 0.1','BLER R 0.1','BER R 0.2','BLER R 0.2','BER R 0.3','BLER R 0.3','BER R 0.4','BLER R 0.4','BER R 0.5','BLER R 0.5','BER R 0.6','BLER R 0.6','BER R 0.7','BLER R 0.7','BER R 0.8','BLER R 0.8','BER R 0.9','BLER R 0.9');
savefig('results_rate_curves/BER_LDPC_960.fig');
saveas(f,'results_rate_curves/BER_LDPC_960.png');

