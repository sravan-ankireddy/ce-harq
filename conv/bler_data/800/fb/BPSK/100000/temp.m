cd all; 

err_thr_grid = 0.0:0.005:0.095;
gs_size = length(err_thr_grid);

SNRdB_low = -10;
SNRdB_high = -4;
SNRdB_step = 0.2;
SNRdB_vec = SNRdB_low:SNRdB_step:SNRdB_high;

snr_data = SNRdB_vec;
bler_data = zeros(gs_size,10,length(snr_data));
ber_data = zeros(gs_size,10,length(snr_data));
ar_data = zeros(gs_size,length(snr_data));

for i = 1:gs_size
    filename = sprintf('fb_data_Conv_800_rate_0.500_err_thr_%.3f_max_rounds_10',err_thr_grid(i));

    data = load(filename);

    bler_data(i,:,:) = data.bler_data;
    ber_data(i,:,:) = data.ber_data;
    ar_data(i,:) = data.ar_data;
end

% Store data 
data_file_name_gs = "fb_data_Conv_800_rate_0.500_err_thr_0.000_to_0.095_max_rounds_10.mat";
save(data_file_name_gs,'ber_data','bler_data','ar_data','snr_data','err_thr_grid');


cd ..;