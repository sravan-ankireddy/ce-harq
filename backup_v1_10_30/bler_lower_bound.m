function eps = bler_lower_bound(N, R, snr_db, channel, modulation, bound)
    % Example : 
    % N = 960; %Blocklength
    % R = 0.9; %Rate
    % snr_db = -5;
    % channel = 'BIAWGN';
    % modulation = 'QPSK';
    % bound = 'conv';
    % bler = bler_lower_bound(N, R, snr_db, channel, modulation, bound);

    % first get rough estimate
    eps = double(5);
    rate = 1.01;
    while rate > R && eps>0
        eps = eps/10;
        if strcmp(bound ,'ach')
            rate = achievable_rate(eps, N, snr_db, channel, modulation);
        elseif strcmp(bound,'conv')
            rate = converse_rate(eps, N, snr_db, channel, modulation);
        end
%         fprintf('eps: %f %d\n',eps, rate);
    end
    % eps is less than BLER
    del_eps = eps/10;
    times = 5;
    while rate<=R
        eps = eps+del_eps;
        if strcmp(bound ,'ach')
            rate1 = achievable_rate(eps, N, snr_db, channel, modulation);
        elseif strcmp(bound ,'conv')
            rate1 = converse_rate(eps, N, snr_db, channel, modulation);
        end

        if rate1<R
            rate = rate1;
        else
            if times>=0
                eps = eps - del_eps;
                del_eps = del_eps/10;
                times = times-1;
%                 fprintf('Once done\n, eps= %f', eps);
            else
                break;
            end

        end 
%         fprintf('eps_2: %f\n, rate : %f',eps, rate);

    end
    eps = eps - del_eps;
end

