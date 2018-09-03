function sweep_results = TG_sweep_1d(var_str, var_vec, total_runs)
%TG_sweep_1d - Preforms 1d sweep over a given variable
%
% Syntax: sweep_results = myFun(input)
%
% Long description



%% Default parameter values
N = 400;  % Total number of neurons
f = 0.5;  % Proportion of excitatory nodes 
tau = 1/4;  % Decay rate of system

sigmaev = 5;  % excitatory standard deviation
sigmaiv = 1;  % inhibitory standard deviation

mue = 0;  % excitatory mean
mui = -mue;  % inhibitory mean

frac_EV_TG = 0.5;  % fraction of eigenvalues used for TG

t_min = 0;
t_max = 10;
t_step = 0.1;



%% Preallocate data array
sweep_results = zeros(length(var_vec), total_runs, 3);


%% Run solvers
% loop over sweep values
for i = 1:length(var_vec)

    % Load input sweep parameter
    eval([var_str ' = ' num2str(var_vec(i)) ';']);
    disp([var_str ' = ' num2str(var_vec(i))]); 

    % Set up parameters
    sigma_vec = [sigmaev, sigmaiv];
    mu_vec = [mue, mui];
    t_vec = t_min:t_step:t_max;


    % loop over runs
    for j = 1:total_runs

        disp(['  run = ' num2str(j)])
        seed = 100 * j;
        
        % run solver
        [V, omega, ~] = ...
            TG_get_eig_matrix(N, f, tau, sigma_vec, mu_vec, seed);

        [~, G_stats] = TG_get_max_growth(t_vec, V, omega, frac_EV_TG);

        % extract data
        sweep_results(i, j, 1) = G_stats.G_max;
        sweep_results(i, j, 2) = G_stats.t_opt;
        sweep_results(i, j, 3) = G_stats.G_init_slope;

        clear V omega G_stats;

    end

end