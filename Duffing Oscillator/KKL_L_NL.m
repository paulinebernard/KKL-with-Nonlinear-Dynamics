close all
clear all

Lambda = [-5;-2];
%Lambda = [-6;-4;-2];
dim_z = length(Lambda); 
nbr_zdot = 3;
dim_z_tot = dim_z*nbr_zdot;

dim_x = 2;

names_obs = {'Fast linear KKL'; 'Slow linear KKL'; 'Nonlinear KKL'};

%% Learning T

t_learning = 20/min(abs(real(Lambda)));

x_lim = [-2, 2; -2, 2]; % dim = dim_x * 2

nb_x_gen = 200;

parfor k = 1:dim_x
    x_gen{k} = linspace(x_lim(k,1), x_lim(k,2),nb_x_gen);  % grid of initial conditions
end

% À modifer selon dim_x :
[X1, X2] = ndgrid(x_gen{1}, x_gen{2});
X0 = [X1(:), X2(:)];

z0 = zeros(dim_z_tot,1);

x_data = zeros(length(X0),dim_x);
z_data = zeros(length(X0),dim_z_tot);

parfor k = 1:length(X0)
    [time, X] = ode45(@(t,X) [compute_xdot(X(1:dim_x));compute_zdot(X(dim_x+1:end),h(X(1:dim_x)),Lambda)],[0,t_learning],[X0(k,:)';z0]);
    x_data(k,:) = X(end,1:dim_x);
    z_data(k,:) = X(end,dim_x+1:end);
end


%% Plot Learning T

% figure
% scatter3(x_data(:,1),x_data(:,2),z_data(:,1))
% xlabel('$x_1$','Interpreter','latex')
% ylabel('$x_2$','Interpreter','latex')
% title('$T_1(x)$','Interpreter','latex')
% 
% figure
% scatter3(x_data(:,1),x_data(:,2),z_data(:,2))
% xlabel('$x_1$','Interpreter','latex')
% ylabel('$x_2$','Interpreter','latex')
% title('$T_2(x)$','Interpreter','latex')

%% convergence time and variance estimation

nb_CI_x = 100;
radius = 1;
angles = rand(nb_CI_x,1)*2*pi;

X0 = radius*[cos(angles),sin(angles)];
%load("x0.mat") % for nb_CI_x = 100 and et dim_x = 2

tf = t_learning;

% Compute the maximum error starting from z0 = T(x0) for each initial 
% condition and each observer

tolerance_x = -1;
tolerance_z = -1;

Z0 = zeros(nb_CI_x, dim_z_tot);

parfor indx0 = 1:nb_CI_x

    [~,ind_min] = min(sum(abs(x_data-X0(indx0, :)),2));
    z0 = z_data(ind_min,:);

    Z0(indx0, :) = z0;

    [time, X] = ode45(@(t,X) [compute_xdot(X(1:dim_x));compute_zdot(X(dim_x+1:end),h(X(1:dim_x)),Lambda)],[0,t_learning],[X0(indx0,:)';z0']);
    
    [~, ~, Err_norm_x, Err_norm_z, ~, ~] = compute_xhat_Tx(time, X, x_data, z_data, dim_x, nbr_zdot, dim_z, -1, -1);

    for obs = 1:nbr_zdot
        tolerance_x = max(tolerance_x, max(Err_norm_x{obs}));
        tolerance_z = max(tolerance_z, max(Err_norm_z{obs}));
    end

end

% Now, we have tolerance_x and tolerance_z, that are the thresholds to
% compute the convergence time

margin_x = 1.1;
margin_z = 1.1;

tolerance_x = margin_x * tolerance_x;
tolerance_z = margin_z * tolerance_z;

% Now, compute the trajectories starting from initial conditions with
% initial error = err0

norm_err_0 = 100;
err0 = rand(1, dim_z);
err0 = norm_err_0*err0/norm(err0);
%load("err0.mat") % for dim_z = 2
err0 = kron(ones(1,nbr_zdot),err0);

t_conv_all = zeros(nbr_zdot, 2, nb_CI_x);
variance_all_noise = zeros(nbr_zdot, 2, nb_CI_x);

noise_level = 0.1;
noise_freq = 10;
noise = @(t) noise_level*sin(noise_freq*t);

parfor indx0 = 1:nb_CI_x

    % Trajectories without noise but with initial error, to estimate 
    % convergence time
    
    z0 = Z0(indx0, :) + err0;

    [time, X] = ode45(@(t,X) [compute_xdot(X(1:dim_x));compute_zdot(X(dim_x+1:end),h(X(1:dim_x)),Lambda)],[0,t_learning],[X0(indx0,:)';z0']);
    
    [xhat_array, Tx_array, Err_norm_x, Err_norm_z, t_conv, ~] = compute_xhat_Tx(time, X, x_data, z_data, dim_x, nbr_zdot, dim_z, tolerance_x, tolerance_z);

    Times{indx0} = time;
    X_all{indx0} = X;
    xhat_all{indx0} = cell2mat(xhat_array);
    Tx_all{indx0} = cell2mat(Tx_array);
    Err_norm_x_all{indx0} = cell2mat(Err_norm_x);
    Err_norm_z_all{indx0} = cell2mat(Err_norm_z);
    t_conv_all(:, :, indx0) = cell2mat(t_conv');

    % Trajectories with noise but without initial error, to estimate 
    % error in asymptotic regime

    z0 = Z0(indx0, :);

    [time, X] = ode45(@(t,X) [compute_xdot(X(1:dim_x));compute_zdot(X(dim_x+1:end),h(X(1:dim_x))+noise(t),Lambda)],[0,t_learning],[X0(indx0,:)';z0']);
    
    [xhat_array, Tx_array, Err_norm_x, Err_norm_z, ~, variance_noise] = compute_xhat_Tx(time, X, x_data, z_data, dim_x, nbr_zdot, dim_z, tolerance_x, tolerance_z);

    Times_noise{indx0} = time;
    X_all_noise{indx0} = X;
    xhat_all_noise{indx0} = cell2mat(xhat_array);
    Tx_all_noise{indx0} = cell2mat(Tx_array);
    Err_norm_x_all_noise{indx0} = cell2mat(Err_norm_x);
    Err_norm_z_all_noise{indx0} = cell2mat(Err_norm_z);
    variance_all_noise(:, :, indx0) = cell2mat(variance_noise');

end

%% Display and plots

min_t_conv = min(t_conv_all, [], 3);
max_t_conv = max(t_conv_all, [], 3);
average_t_conv = sum(t_conv_all, 3)./nb_CI_x;

min_variance_noise = min(variance_all_noise, [], 3);
max_variance_noise = max(variance_all_noise, [], 3);
average_variance_noise = sum(variance_all_noise, 3)./nb_CI_x;

for obs = 1:nbr_zdot
    disp(names_obs{obs})
    disp(['Average convergence time : ', num2str(average_t_conv(obs,1)), ' in x  and ', num2str(average_t_conv(obs,2)), ' in z'])
    disp(['Average variance with noise : ', num2str(average_variance_noise(obs,1)), ' in x  and ', num2str(average_variance_noise(obs,2)), ' in z'])
end

close all


% Err_norm = sqrt(sum((X_all{1}(:, 1:dim_x) - xhat_all{1}(:, 1:dim_x)).^2, 2));
% Err_norm_sort = smooth(sort(Err_norm));
% [m, i] = max(diff(Err_norm_sort));
% tolerance_x = Err_norm(i)
% 
% figure;
% plot(Times{1}, Err_norm)
% hold on
% plot([Times{1}(1),Times{1}(end)], [1,1]*tolerance_x)
% figure;
% plot(1:length(Err_norm), Err_norm_sort);
% hold on
% plot([1,length(Err_norm)], [1,1]*tolerance_x)
% figure;
% plot(1:length(Err_norm)-1, diff(Err_norm_sort));

% return
% close all
% figure
% plot(Times{1}, X_all{1}(:,1))
% hold on
% plot(Times{1}, xhat_all{1}(:,1))
% 
% figure
% plot(Times{1}, X_all{1}(:,3:5))
% hold on
% plot(Times{1}, Tx_all{1}(:,1:3))
%

for obs = 1:nbr_zdot
    start = (obs-1)*dim_z+1;
    stop = obs*dim_z;
    start_x = (obs-1)*2+1;
    stop_x = obs*2;
    figure(1)
    plot(Times{1}, Err_norm_x_all{1}(:, obs))
    hold on
    figure(2)
    plot(Times{1}, Err_norm_z_all{1}(:, obs))
    hold on
    % obs
    % t_conv_all(obs, :, 1)
    % ind_conv_all(obs, :, 1)
    % variance_noise_all(obs, :, 1)
end

figure(1)
leg = legend(names_obs, 'interpreter', 'latex');
xlabel('$t$', 'interpreter', 'latex')
%ylabel('$\|x(t)-\hat x(t)\|$', 'interpreter', 'latex')
xlim([0, tf])
myPrintPDF(1,leg,'error_x_conv')

figure(2)
leg = legend(names_obs, 'interpreter', 'latex');
xlabel('$t$', 'interpreter', 'latex')
%ylabel('$\|z(t)-T(x(t))\|$', 'interpreter', 'latex')
xlim([0, tf])
myPrintPDF(2,leg,'error_z_conv')



for obs = 1:nbr_zdot
    start = (obs-1)*dim_z+1;
    stop = obs*dim_z;
    start_x = (obs-1)*2+1;
    stop_x = obs*2;
    figure(3)
    plot(Times_noise{1}, Err_norm_x_all_noise{1}(:, obs))
    hold on
    figure(4)
    plot(Times_noise{1}, Err_norm_z_all_noise{1}(:, obs))
    hold on
    % obs
    % t_conv_all(obs, :, 1)
    % ind_conv_all(obs, :, 1)
    % variance_noise_all(obs, :, 1)
end

figure(3)
leg = legend(names_obs, 'interpreter', 'latex');
xlabel('$t$', 'interpreter', 'latex')
%ylabel('$\|x(t)-\hat x(t)\|$', 'interpreter', 'latex')
xlim([0, tf])
myPrintPDF(3,leg,'error_x_noise')

figure(4)
leg = legend(names_obs, 'interpreter', 'latex');
xlabel('$t$', 'interpreter', 'latex')
%ylabel('$\|z(t)-T(x(t))\|$', 'interpreter', 'latex')
xlim([0, tf])
myPrintPDF(4,leg,'error_z_noise')

colors = ["#000000";"#0072BD";"#D95319";"#EDB120"];
figure(5)
plot(Times{1}, X_all{1}(:, 1:2),'Color',colors{1},'Linewidth',2)
hold on
h = zeros(nbr_zdot+1, 1);
h(1) = plot(NaN,NaN,'Color',colors{1});
for obs = 1:nbr_zdot
    start_x = (obs-1)*2+1;
    stop_x = obs*2;
    plot(Times{1}, xhat_all{1}(:, start_x:stop_x),'Color',colors{obs+1},'Linewidth',2)
    h(obs+1) = plot(NaN,NaN,'Color',colors{obs+1});
    hold on
end
leg = legend(h, ['System';names_obs],'Fontsize',9);

