clear; close all; clc;
noTransportBlocks = 10e2;
SNRdB_vec = -10:0.2:0; % SNR in dB
BLER = zeros(2,length(SNRdB_vec));
noNewBlocks = zeros(2,length(SNRdB_vec));
targetCodeRate = 0.9;
rvSeq_IR = [0 2 3 1 0 2 3 1 0 2];
rvSeq_Chase = zeros(size(rvSeq_IR));
res_folder = "ref_curves";
tic;
parpool_size = min(64,feature('numcores'));
if (parpool_size ~= 8 && parpool_size ~= 10)
    parpool(parpool_size);
end
rvSeq = rvSeq_IR;
codeLen = -1*zeros(size(SNRdB_vec));
codeRate = -1*zeros(size(SNRdB_vec));
for j = 1:2
    if (j == 1) % Chase
        rvSeq = rvSeq_Chase;
    else
        rvSeq = rvSeq_IR;
    end
    
    parfor i = 1:length(SNRdB_vec)
        SNRdB = SNRdB_vec(i);
    
        rng("default");
        
        % DL-SCH parameters
        NHARQProcesses = 1; % Number of parallel HARQ processes to use
            
        % Create DL-SCH encoder object
        encodeDLSCH = nrDLSCH;
        encodeDLSCH.MultipleHARQProcesses = true;
        encodeDLSCH.TargetCodeRate = targetCodeRate;
        
        % Create DL-SCH decoder object
        decodeDLSCH = nrDLSCHDecoder;
        decodeDLSCH.MultipleHARQProcesses = true;
        decodeDLSCH.TargetCodeRate = targetCodeRate;
        decodeDLSCH.LDPCDecodingAlgorithm = "Normalized min-sum";
        decodeDLSCH.MaximumLDPCIterationCount = 6;
        
        % Numerology
        SCS = 15;                         % SCS: 15, 30, 60, 120 or 240 (kHz)
        NRB = 3;                         % BW in number of RBs (52 RBs at 15 kHz SCS for 10 MHz BW)
        
        carrier = nrCarrierConfig;
        carrier.NSizeGrid = NRB;
        carrier.SubcarrierSpacing = SCS;
        carrier.CyclicPrefix = "Normal";  % "Normal" or "Extended"
        
        modulation = "QPSK";             % Modulation scheme
        
        pdsch = nrPDSCHConfig;
        pdsch.Modulation = modulation;
        pdsch.PRBSet = 0:NRB-1;           % Assume full band allocation
        pdsch.NumLayers = 1;              % Assume only one layer and one codeword
        
        harqEntity = HARQEntity(0:NHARQProcesses-1,rvSeq,pdsch.NumCodewords);
        
        % Initialize loop variables
        noiseVar = 1./(10.^(SNRdB/10)); % Noise variance
        numBlkErr = 0;                  % Number of block errors
        numRxBits = [];                 % Number of successfully received bits per transmission
        txedTrBlkSizes = [];            % Number of transmitted info bits per transmission
        trBlk = 0;
        for nTrBlk = 1:noTransportBlocks
            % A transport block or transmission time interval (TTI) corresponds to
            % one slot
            carrier.NSlot = carrier.NSlot+1;
        
             % Generate PDSCH indices info, which is used to calculate the transport
            % block size
            [~,pdschInfo] = nrPDSCHIndices(carrier,pdsch);
        
            % Calculate transport block sizes
            Xoh_PDSCH = 0;
            trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet),pdschInfo.NREPerPRB,targetCodeRate,Xoh_PDSCH);
        
            % Get new transport blocks and flush decoder soft buffer, as required
            for cwIdx = 1:pdsch.NumCodewords
                if harqEntity.NewData(cwIdx)
                    % Create and store a new transport block for transmission
                    trBlk = randi([0 1],trBlkSizes(cwIdx),1);
                    setTransportBlock(encodeDLSCH,trBlk,cwIdx-1,harqEntity.HARQProcessID);
        
                    % If the previous RV sequence ends without successful decoding,
                    % flush the soft buffer explicitly
                    if harqEntity.SequenceTimeout(cwIdx)
                        resetSoftBuffer(decodeDLSCH,cwIdx-1,harqEntity.HARQProcessID);
                    end
                end
            end
        
            codedTrBlock = encodeDLSCH(pdsch.Modulation,pdsch.NumLayers,pdschInfo.G, ...
                harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
                
            modOut = nrPDSCH(carrier,pdsch,codedTrBlock);
        
            rxSig = awgn(modOut,SNRdB);
        
            rxLLR = nrPDSCHDecode(carrier,pdsch,rxSig,noiseVar);
            
            codeLen(i) = length(codedTrBlock);
            codeRate(i) =length(trBlk)/codeLen(i);
            decodeDLSCH.TransportBlockLength = trBlkSizes;
            [decbits,blkerr] = decodeDLSCH(rxLLR,pdsch.Modulation,pdsch.NumLayers, ...
                harqEntity.RedundancyVersion,harqEntity.HARQProcessID);

            % Store values to calculate throughput (only for active transport blocks)
            if(any(trBlkSizes ~= 0))
                numRxBits = [numRxBits trBlkSizes.*(1-blkerr)];
                txedTrBlkSizes = [txedTrBlkSizes trBlkSizes];
            end
    
            statusReport = updateAndAdvance(harqEntity,blkerr,trBlkSizes,pdschInfo.G);
    
            if blkerr && harqEntity.NewData
                numBlkErr = numBlkErr + 1;
            end
    
            if harqEntity.NewData
                noNewBlocks(j,i) = noNewBlocks(j,i) + 1;
            end
        
            disp("Slot "+(nTrBlk)+". "+statusReport);
           
        end % for nTrBlk = 1:noTransportBlocks
        
        totalNumRxBits = sum(numRxBits,2);   % Number of successfully received bits
        
        disp("Block Error Rate: "+string(numBlkErr/noNewBlocks(j,i)))
        
        BLER(j,i) = numBlkErr/noNewBlocks(j,i);
    end
end
toc;

actualCodeRate = codeRate(end);
%%
figure('Renderer','painters','Position',[1000 400 800 500]);
f = semilogy(SNRdB_vec,squeeze(BLER(1,:)),'b-o');
hold on;
semilogy(SNRdB_vec,squeeze(BLER(2,:)),'r-d');
fs = 16;
xlabel('SNR','FontSize',fs);
ylabel('BLER','FontSize',fs);
max_rounds = length(rvSeq);
leg_IR = sprintf('HARQ-IR Rate %.3f, max. %d rounds',actualCodeRate, max_rounds);
leg_Chase = sprintf('HARQ-Chase Rate %.3f, max. %d rounds',actualCodeRate, max_rounds);
legend(leg_Chase,leg_IR,'Location','southwest','FontSize',fs);
title_name = sprintf('HARQ IR vs Chase combining : LDPC %d, Rate %.3f, max. %d rounds',codeLen(end), actualCodeRate, max_rounds);
title(title_name,'FontSize',fs);
if ~exist(res_folder,'dir')
    mkdir(res_folder);
end
fig_name = sprintf('ref_curves/IRvsChase_BLER_LDPC_%d_HARQ_%d_numF_%d.fig',codeLen(end),max_rounds,noTransportBlocks);
savefig(fig_name);
png_name = sprintf('ref_curves/IRvsChase_BLER_LDPC_%d_HARQ_%d_numF_%d.png',codeLen(end),max_rounds,noTransportBlocks);
saveas(f,png_name);

if (parpool_size ~= 8 || parpool_size ~= 10)
    delete(gcp('nocreate'));
end