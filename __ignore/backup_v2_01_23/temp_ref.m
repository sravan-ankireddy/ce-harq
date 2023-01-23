function out = temp_ref(channel,SNRdB,modulation,N,K,R,PHY_code,MAC_code,feedback_mode,combining_scheme,rvSeq,ncb,Nref,max_iter,nlayers,dec_type,data,rxLLR,data_est,err_thr,err_thr_ada_list_est,err_thr_ada_scheme,i_s,max_rounds,counts,num_err,comm_mod,mod_approx,seed)
    
    rng(seed);
    int_state = seed;
    base_rate = R;
    k = bits_per_symbol(modulation);
    M = 2^k;
    Avg_rounds_FB = 0;
    rxLLR_HARQ_buffer = rxLLR;
    rxLLR_FB_buffer = [];
    rxLLR_FB_buffer_outer = [];
    data_est_FB = data_est;
    data_est_FB_prev = data_est;
    decision_switch = 0;
    fb_scheme = "HARQ";
    num_err_FB_prev = num_err;
    num_err_vec = zeros(max_rounds,1);
    num_err_vec(1) = num_err;

    noiseVar = 1./(10.^(SNRdB/10));
    rep_coding = 0;
    tc_rp = 0.2;

    comp_rates = 1./(2:12);
    minR = comp_rates(end);
    rc = 0;

    % FIX ME
    rv = 0;
    rvSeq_ind = 1;

    NREPerPRB = 12*4; % For URLLC, 2-7 is the typical choice

    min_bler = 1e-4;
    harq_data_path = sprintf('bler_data/awgn/%d/%s/harq/%s/100000/harq_data_Conv_%d_rate_0.833_rate_0.083_max_rounds_%d.mat', N, dec_type, modulation, N, max_rounds);
    acomp_table_path = sprintf('lut_data/acomp_%d_ns_100000.mat', K);
    if (err_thr_ada_scheme == "est")
        harq_data = load(harq_data_path);
        acomp_table = load(acom_table_path);
    end
    prev_reset_round = 0;%max_rounds;

    
    for i_r = 1:max_rounds-1

        num_err_FB = sum(data ~= double(data_est_FB));
        err_per = num_err_FB/K;

        % Adaptive err thr
        if (err_thr_ada_scheme == "est")
            if (fb_scheme == "FB")
                err_thr = err_thr_select(harq_data,acomp_table,targetErrCodeRate,SNRdB,max_rounds - prev_reset_round,max_rounds - i_r,min_bler);
            else
                err_thr = err_thr_select(harq_data,acomp_table,R,SNRdB,max_rounds,max_rounds - i_r,min_bler);
            end
        end

        % Once the error becomes sparse enough, stay on FB scheme
        % Within the FB scheme, you can choose FB/HARQ based on error
        % vector
        if (decision_switch == 0)
            if (err_per <= err_thr && err_thr > 0)
                % disp("found sparse err ")
                % Also check if the error is compressible
                data_est_err_temp = mod(data+double(data_est_FB),2);
                % iterative counts
                p_cur = max(1,round(100 * sum(data_est_err_temp)/length(data_est_err_temp)));
                counts_temp = [100-p_cur p_cur];
                err_seq_temp = arithenco(data_est_err_temp+1,counts_temp);
                % err_seq_temp = arithenco(data_est_err_temp+1,counts);
                if (length(err_seq_temp) < length(data_est_err_temp))
                    fb_scheme = "FB";
                    decision_switch = 1;
                    counts = counts_temp;
                    % disp(" switched to feedback ")
                else
                    fb_scheme = "HARQ";
                end
            else
                fb_scheme = "HARQ";
            end
        end

        % Update the counter
        Avg_rounds_FB = Avg_rounds_FB + 1; 
        
        % Stay on HARQ if error is not sparse yet
        if (fb_scheme == "FB")
            % Compress the error assuming free feedback
            data_est_err = mod(data+double(data_est_FB),2);

            % Choose HARQ/FB based on error vector
            num_err_FB_cur = sum(data_est_err);
            if (num_err_FB_cur > num_err_FB_prev)
                % No. of errors increased from prev. round; use HARQ on FB
                inner_scheme = "HARQ";

                % Reset latest est to prev round est
                data_est_FB = data_est_FB_prev;
                num_err_FB = num_err_FB_prev;
            else
                inner_scheme = "FB";
                % Empty the chase buffer if prev round reduced the errors

                rxLLR_FB_buffer = [];
                rxLLR_FB_buffer_outer = [];
                % Reset rvSeq start
                rvSeq_ind = 1;
            end

            %% Tx
            if (inner_scheme == "FB")                
                % iterative counts
                p_cur = max(1,round(100 * sum(data_est_err)/length(data_est_err)));
                counts = [100-p_cur p_cur];
                err_seq = arithenco(data_est_err+1,counts);
    
                assert(length(err_seq) < length(data_est_err),'Compression failed in FB round %d',i_r);

                % PHY-MAC mode
                if (feedback_mode == "MAC")

                    % MAC layer coding -- inner code : compress the error and perform FEC
                    if (MAC_code == "Conv")
                        % find the smallest rate for compression : coding will be only for comp -> K, not N
                        targetErrCodeRate = (length(err_seq)+3)/K;
                        R_mac = comp_rate_conv(targetErrCodeRate, minR, comp_rates);

                        % encode to length <= K
                        comp_code_inner = conv_enc(err_seq, R_mac);

                        % interleaver between inner and outer code
                        comp_code_inner = randintrlv(comp_code_inner,int_state);

                        % append zeros to match length to K
                        nz = K - length(comp_code_inner);
                        comp_code_inner = [comp_code_inner; zeros(nz,1)];  
						
						% standardizing notation
						K_mac = K;
						N_mac = N;
                    end
				
					% PHY layer coding -- perform outer coding with rate K/N
					if (PHY_code == "Conv")
                        comp_code_outer = conv_enc(comp_code_inner, base_rate);
					end
				
				% Direct PHY+MAC coding with lowest possible rate : PHY_code is the only code that matters
				else
					if (PHY_code == "Conv")
						% find the smallest rate for compression : coding will be only for comp -> K, not N
						targetErrCodeRate = (length(err_seq)+3)/N;
						R_phy = comp_rate_conv(targetErrCodeRate, minR, comp_rates);

						% encode to length <= K
						comp_code_inner = conv_enc(err_seq, R_phy);

						% interleaver between inner and outer code
						comp_code_inner = randintrlv(comp_code_inner,int_state);

						comp_code_outer = comp_code_inner;
					end
                end % ending PHY-MAC if 
                
				dataIn = comp_code_outer;
            end % Stay on HARQ if inner_scheme is not FB; no new dataIn generated

            % generate llrs
            % Pass through channel, modulation, democulation
            rxLLR_FB = transmit_data (channel, dataIn, SNRdB, modulation, dec_type);

            % Chase Combining
            if (dec_type == "hard")
                rxLLR_FB_buffer = [rxLLR_FB_buffer rxLLR_FB];
                rxLLR_FB = round(mean(rxLLR_FB_buffer, 2));
            else
                rxLLR_FB_buffer = [rxLLR_FB_buffer rxLLR_FB];
                rxLLR_FB = sum(rxLLR_FB_buffer, 2);
            end

            % Decoding
            if (PHY_code == "Conv")
                if (feedback_mode == "MAC")
                    % first decode outer code
                    outer_err_seq_est = conv_dec(rxLLR_FB,base_rate,dec_type);

                    % remove the zero padding
                    outer_err_seq_est = outer_err_seq_est(1:end-nz);
                    
                    % deinterleaver between inner and outer code
                    outer_err_seq_est = randdeintrlv(outer_err_seq_est,int_state);
                    
                    % next decode inner : map to 1 - 2*c for soft
                    if (dec_type == "hard")
                        rxLLR_FB = outer_err_seq_est;
                    else
                        rxLLR_FB = 1 - 2*outer_err_seq_est;
                    end
                    inner_err_seq_est = conv_dec(rxLLR_FB,comp_rate,dec_type);
                else
                    inner_err_seq_est = conv_dec(rxLLR_FB,comp_rate,dec_type);
                end
            end

            % Decompress
            err_deseq_est = arithdeco(double(inner_err_seq_est),counts,length(data))-1;

            % Correct regardless of CRC but store the data (before
            % correction) : except for last round
            data_est_FB_prev = data_est_FB;
            num_err_FB_prev = sum(data ~= double(data_est_FB));

            %% FIX ME : temp 
            data_est_FB_temp = mod(data_est_FB + err_deseq_est,2);
            num_err_FB_temp = sum(data ~= double(data_est_FB_temp));
            if ((i_r < max_rounds - 1) || (i_r == max_rounds -1 && num_err_FB_temp == 0))
                data_est_FB = mod(data_est_FB + err_deseq_est,2);
                num_err_FB = sum(data ~= double(data_est_FB));
            end
            
            % update err count
            num_err_vec(i_r+1) = num_err_FB;
            
            rvSeq_ind = rvSeq_ind + 1;
            
            % Break if crc passes
            if (num_err_FB == 0)
%                 disp(" ************* CORRECTED WITH SPARSE ************* ")
                break;
            end
        else
            % Retransmit
            % generate new tx data
            
            if (PHY_code == "Conv")         
                [dataIn, ~] = conv_enc(data, R);
            end

            newRxLLR_HARQ = transmit_data (channel, dataIn, SNRdB, modulation, dec_type);

            if (dec_type == "hard")
                rxLLR_HARQ_buffer = [rxLLR_HARQ_buffer newRxLLR_HARQ];
                rxLLR_FB_HARQ = round(mean(rxLLR_HARQ_buffer,2));
            else
                rxLLR_HARQ_buffer = [rxLLR_HARQ_buffer newRxLLR_HARQ];
                rxLLR_FB_HARQ = sum(rxLLR_HARQ_buffer,2);
            end
            
            % Rate recovery and Decoding
            if (PHY_code == "Conv")
                data_est_FB = conv_dec(rxLLR_FB_HARQ, R, dec_type);
            end

            % Check for errors : post descrambling
            num_err_FB = sum(mod(data+double(data_est_FB),2));

            % update err count
            num_err_vec(i_r+1) = num_err_FB;
            
            % Break if crc passes
            if (num_err_FB == 0)
                % disp(" ************* CORRECTED WITH HARQ ************* ")
                break;
            end
        end

    end

    % output
    out.Avg_rounds_FB = Avg_rounds_FB;
    out.num_err_FB = num_err_FB;
    out.num_err_vec = num_err_vec;

    function comp_rate = comp_rate_conv(targetErrCodeRate, minR, comp_rates)

        if (targetErrCodeRate >= minR)
            % get the smallest rate possible
            [~, min_ind] = min(abs(targetErrCodeRate - comp_rates));
            if (comp_rates(min_ind) > targetErrCodeRate)
                comp_rate = comp_rates(min_ind);
            else
                while (comp_rates(min_ind) < targetErrCodeRate)
                    min_ind = min_ind - 1;
                end
                comp_rate = comp_rates(min_ind);
            end
        else
            comp_rate = comp_rates(end);
            min_ind = length(comp_rates) - 1;
        end
    end


    function generate_llr = transmit_data (channel, dataIn, SNRdB, modulation, dec_type)

		% Symbol Modulation
		if (modulation == "BPSK")
			bpskModulator = comm.BPSKModulator;
			bpskModulator.PhaseOffset = pi/4;
			txSig = bpskModulator(dataIn);
		else
			txSig = qammod(dataIn,M,'InputType','bit','UnitAveragePower',true);
		end

		if (channel == "awgn")
			newRxSig = awgn(txSig,SNRdB,'measured');
		elseif (channel == "rayleigh")
			h = sqrt(rand(1)^2 + rand(1)^2);
			newRxSig = awgn(h*txSig,SNRdB,'measured');
		end

		% Symbol demod
		if (modulation == "BPSK")
			if (dec_type == "hard")
				bpskDemodulator = comm.BPSKDemodulator; 
				bpskDemodulator.PhaseOffset = pi/4; 
				bpskDemodulator.DecisionMethod = 'Hard decision';
				generate_llr = bpskDemodulator(newRxSig);
			else
				bpskDemodulator = comm.BPSKDemodulator; 
				bpskDemodulator.PhaseOffset = pi/4; 
				bpskDemodulator.DecisionMethod = 'Approximate log-likelihood ratio';
				generate_llr = bpskDemodulator(newRxSig);
			end
		else
			if (dec_type == "hard")
				generate_llr = qamdemod(newRxSig, M, OutputType='bit', UnitAveragePower=true);
			else
				generate_llr = qamdemod(newRxSig, M, OutputType='approxllr', UnitAveragePower=true, NoiseVariance=noiseVar);
			end
        end
    end

end