function [xhat_array, Tx_array, Err_norm_x, Err_norm_z, t_conv, variance_noise] = compute_xhat_Tx(time, X, x_data, z_data, dim_x, nbr_zdot, dim_z, tolerance_x, tolerance_z)

% X is a trajectory of (x,z): size = (temps de simu) * (dim x + dim_z_tot)

nb_time = length(time);

% Tx = zeros(nb_time,dim_z_tot);
% xhat = zeros(nb_time,dim_x*nbr_zdot);
% t_conv = zeros(nbr_zdot, 2);
% ind_conv = zeros(nbr_zdot, 2);

parfor obs = 1:nbr_zdot

    start = (obs-1)*dim_z+1;
    stop = obs*dim_z;

    xhat_array{obs} = zeros(nb_time, dim_x);
    Tx_array{obs} = zeros(nb_time, dim_z);

    for ind = 1:nb_time

        [~,ind_min] = min(sum(abs(z_data(:,start:stop)-X(ind,(dim_x+start):(dim_x+stop))),2));
        xhat_array{obs}(ind,:) = x_data(ind_min,:);

        [~,ind_min] = min(sum(abs(x_data-X(ind,1:dim_x)),2));
        Tx_array{obs}(ind,:) = z_data(ind_min,start:stop);

    end

    Err_norm_x{obs} = sqrt(sum((X(:,1:dim_x) - xhat_array{obs}).^2, 2));
    Err_norm_z{obs} = sqrt(sum((X(:,(dim_x+start):(dim_x+stop)) - Tx_array{obs}).^2, 2));

    t_conv{obs} = [Inf, Inf];
    variance_noise{obs} = [Inf, Inf];

    if tolerance_x > 0
        t = conv_time(time, X(:, 1:dim_x), xhat_array{obs}, tolerance_x);
        t_conv{obs}(1) = t;
    end
    if tolerance_z > 0
        t = conv_time(time, X(:, (dim_x+start):(dim_x+stop)), Tx_array{obs}, tolerance_z);
        t_conv{obs}(2) = t;
    end
    variance_noise{obs}(1) = max(Err_norm_x{obs})
    variance_noise{obs}(2) = max(Err_norm_z{obs});
end

end