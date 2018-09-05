function TG_main(total_runs, var_str, var_vec)
%TG_main - Main function used to run demos and sweeps
%
% Syntax: TG_main(total_runs, var_str, var_vec)
% 
% Data saved to 'output/TG_main/'
%
% Input: 
%   total_runs = number of draws to do
%     - if total_runs = 0, then run demo and output omega and G plots
%     - otherwise run a sweep over var_str
%   var_str = cell of strings of valiable names to change
%     - if total_runs > 0, only supports 1d sweeps
%   var_vec = vector of values
%     - if total_runs = 0, input vector such that 
%         length(var_vec) == length(var_str)
% 
% >> TG_main  % VALID: runs demo for default values with plots 
% >> TG_main(0)  % VALID: runs demo for default values with plots 
% >> TG_main(5)  % INVALID: sweep requested but no paramter given
% >> TG_main(0, {'N'}, [150])  % VALID
% >> TG_main(0, {'N', 'tau'}, [150, 1/5])  % VALID
% >> TG_main(5, {'N'}, [150, 160, 170])  % VALID



%% check input and init
addpath(genpath('src'))
close all;

if nargin == 0
    total_runs = 0;
    var_str = {};
    var_vec = [];
end

if nargin == 1
    var_str = {};
    var_vec = [];
end

if (total_runs > 0) && (isempty(var_str))
    error('Sweep requested but no paramter given!')
end

%% load default parameters
N = 100;  % Total number of neurons
f = 0.5;  % Proportion of excitatory nodes 
tau = 1/4;  % Decay rate of system

sigmae = 5;  % excitatory standard deviation
sigmai = 1;  % inhibitory standard deviation

mue = 0;  % excitatory mean
mui = -mue;  % inhibitory mean

frac_EV_TG = 0.5;  % fraction of eigenvalues used for TG

t_min = 0;
t_max = 10;
t_step = 0.1;

%% Run solver for default values if no input is given
%% Modify single value, run solver and plot 
if (total_runs == 0) 
    disp('Demo mode, 1 run for default parameters with plots')

    % Set output directory 
    out_dir = 'output/TG_main/demo/';
    
    % Load changed parameters if given
    if ~isempty(var_str)
        disp('Parameters modified')
        for i = 1:length(var_str)
            eval([var_str{i} ' = ' num2str(var_vec(i)) ';']);
        end
    end

    % set input for solver
    seed = 0;
    sigma_vec = [sigmae, sigmai];
    mu_vec = [mue, mui];
    t_vec = t_min:t_step:t_max;

    % Save input to file
    file = [out_dir 'results.txt'];
    fileID = fopen(file, 'w');
    

    var_names = {'N', 'f', 'tau', 'sigmae', 'sigmai', 'mue', 'mui', ...
        'frac_EV_TG', 't_min', 't_max', 't_step', 'seed'};

    fprintf(fileID, datestr(now,'mm/dd/yyyy HH:MM:SS\n'));
    fprintf(fileID, '1 run, default values, plots on\n\n');
    fprintf(fileID, 'Parameter values loaded : \n');
    fprintf(fileID, '---------------------------\n');
    for i = 1:length(var_names)
        name = var_names{i};
        eval(['val = ' name ';']);
        fprintf(fileID, '  %-12s = %g\n',name,val);
    end
    fprintf(fileID, '---------------------------\n\n');

    % Get eigenvalues of J
    [V, omega, J_mat] = ...
        TG_get_eig_matrix(N, f, tau, sigma_vec, mu_vec, seed);

    % Get max growth
    [G_vec, G_stats] = TG_get_max_growth(t_vec, V, omega, frac_EV_TG);

    % Write G stats to file
    fprintf(fileID, 'Results : \n');
    fprintf(fileID, '---------------------------\n');
    fprintf(fileID,'  %-12s = %g\n','G_max', G_stats.G_max);
    fprintf(fileID,'  %-12s = %g\n','t_opt', G_stats.t_opt);
    fprintf(fileID,'  %-12s = %g\n','G_init_slope', G_stats.G_init_slope);
    
    % Plot eigenvalues
    TG_plot_eigvals(N, f, tau, sigma_vec, mu_vec, omega);
    plot_export_fig(-1, [out_dir 'figures/plot_eigvals'], 14, 7/5, 18);

    % Plot max growth
    TG_plot_max_growth(G_vec, t_vec, G_stats)
    plot_export_fig(-1, [out_dir 'figures/plot_max_growth'], 14, 7/5, 18);

    % Save data to dir
    save([out_dir 'results.mat']);


else 
    disp('Sweep ')

    % Check input for 1d sweep
    if length(var_str) > 1
        error('Only supports 1d sweeps')
    end

    % Convert single cell to str
    var_str = var_str{1};

    % Set output directory 
    out_dir = 'output/TG_main/1dsweep/';
    
    % Save input to file
    file = [out_dir 'results.txt'];
    fileID = fopen(file, 'w');

    var_names = {'N', 'f', 'tau', 'sigmae', 'sigmai', 'mue', 'mui', ...
        'frac_EV_TG', 't_min', 't_max', 't_step'};

    fprintf(fileID, datestr(now,'mm/dd/yyyy HH:MM:SS\n'));

    fprintf(fileID,'%g runs, sweep over %s, plots off\n\n',total_runs,var_str);

    fprintf(fileID, 'Default parameter values loaded : \n');
    fprintf(fileID, '---------------------------\n');
    for i = 1:length(var_names)
        name = var_names{i};
        eval(['val = ' name ';']);
        fprintf(fileID, '  %-12s = %g\n',name,val);
    end
    fprintf(fileID, '---------------------------\n\n');

    fprintf(fileID, 'Sweep parameter values loaded : \n');
    fprintf(fileID, '---------------------------\n');
    fprintf(fileID, '    %s \n', var_str);
    for i = 1:length(var_vec)
        fprintf(fileID, '    %g \n', var_vec(i));
    end
    fprintf(fileID, '---------------------------\n\n');


    % Preallocate data array
    sweep_results = zeros(length(var_vec), total_runs, 3);

    %% Run solvers
    % loop over sweep values
    for i = 1:length(var_vec)

        % Load input sweep parameter
        eval([var_str ' = ' num2str(var_vec(i)) ';']);
        disp([var_str ' = ' num2str(var_vec(i))]); 

        % Set up parameters
        sigma_vec = [sigmae, sigmai];
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

            % Write G stats to file
            fprintf(fileID, 'Results (%s = %g, %g): \n', var_str, var_vec(i), j);
            fprintf(fileID, '---------------------------\n');
            fprintf(fileID,'  %-12s = %g\n','G_max', G_stats.G_max);
            fprintf(fileID,'  %-12s = %g\n','t_opt', G_stats.t_opt);
            fprintf(fileID,'  %-12s = %g\n','G_init_slope', G_stats.G_init_slope);
            fprintf(fileID, '---------------------------\n\n');
            
            % extract data
            sweep_results(i, j, 1) = G_stats.G_max;
            sweep_results(i, j, 2) = G_stats.t_opt;
            sweep_results(i, j, 3) = G_stats.G_init_slope;

            clear V omega G_stats;

        end
    end

    % TODO: plots??


    % Save data to dir
    save([out_dir 'results.mat']);

end


disp('Done!')

end





