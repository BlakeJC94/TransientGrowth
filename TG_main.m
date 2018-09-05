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

if ischar(var_str)
    var_str = {var_str};
end

if (total_runs > 0) && (isempty(var_str))
    error('Sweep requested but no paramter given!')
end



%% load default parameters
Params = TG_main_parameters;
display(Params);

%% Run solver for default values if no input is given
%% Modify single value, run solver and plot 
if (total_runs == 0) 
    disp('Demo mode, 1 run for default parameters with plots')

    % Set output directory 
    out_dir = 'output/TG_main/demo/';
    
    % Load changed parameters if given
    if ~isempty(var_str)

        % Throw error if multiple values given for var_str
        if length(var_str) ~= length(var_vec)
            error('Parameter name value mismatch in demo mode!')
        end

        disp('Parameters modified')
        for i = 1:length(var_str)
            eval(['Params.' var_str{i} ' = ' num2str(var_vec(i)) ';']);
        end
    end
    
    % Save input to file
    TG_write_input(Params, out_dir);

    % Get eigenvalues of J
    [V, omega, J_mat] = TG_get_eig_matrix(Params);  %#ok<ASGLU>

    % Get max growth
    [G_vec, G_stats] = TG_get_max_growth(Params, V, omega);

    % Write G_stats to file
    TG_write_output(G_stats, out_dir);
    
    % Plot eigenvalues
    TG_plot_eigvals(Params, omega, out_dir);

    % Plot max growth
    TG_plot_max_growth(Params, G_vec, G_stats, out_dir);

    % Save data to dir
    save([out_dir 'results.mat']);



else 

    disp('Sweep mode')

    % Check input for 1d sweep
    if length(var_str) > 1
        error('Only supports 1d sweeps')
    end

    % Convert single cell to str
    var_str = var_str{1};

    % Set output directory 
    out_dir = 'output/TG_main/1dsweep/';
    
    % Save input to file
    TG_write_input(Params, out_dir);

    % Preallocate data array
    sweep_results = zeros(length(var_vec), total_runs, 3);

    %% Run solvers
    % loop over sweep values
    disp('=== Begin solving ===')
    for i = 1:length(var_vec)

        % Load input sweep parameter
        eval(['Params.' var_str ' = ' num2str(var_vec(i)) ';']);
        disp([var_str ' = ' num2str(var_vec(i))]); 

        % loop over runs
        for j = 1:total_runs

            disp(['  run = ' num2str(j)])
            Params.seed = 100 * j;
            
            % run solver
            [V, omega, ~] = TG_get_eig_matrix(Params);
            [~, G_stats] = TG_get_max_growth(Params, V, omega);

            % Write G stats to file
            TG_write_output(G_stats, out_dir, var_str, var_vec(i), j);
            
            % extract data
            sweep_results(i, j, 1) = G_stats.G_max;
            sweep_results(i, j, 2) = G_stats.t_opt;
            sweep_results(i, j, 3) = G_stats.G_init_slope;

            clear V omega G_stats;

        end
    end
    disp('=== End solving ===')

    % Save data to dir
    save([out_dir 'results.mat']);

    % Plot results 
    TG_plot_1d_sweep(sweep_results, var_str, var_vec, out_dir)

end


disp('Done!')

end





