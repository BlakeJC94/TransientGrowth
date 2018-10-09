function TG_main_sweep_1d(total_runs, var_str, var_vec)
%TG_main - Main function used to run demos and sweeps
%
% Syntax: TG_main_sweep_1d(total_runs, var_str, var_vec)
% 
% e.g. >> TG_main_sweep_1d(5, 'N', [150, 160, 170])  
%
% Input: 
%   total_runs : number of draws to do
%   var_str : string of variable name to change
%   var_vec = vector of values
% 
% Data saved to 'output/TG_main_sweep_1d/'


%% Initialise

addpath(genpath('src'))
close all;

% check input 
if nargin < 3
    error('Requires 3 input arguments')
end

% Set output directory 
out_dir = 'output/TG_main_sweep_1d/';




%% Load parameters

% load default parameters
Params = TG_parameters;
display(Params);

% Save default input to file
TG_write_input(Params, out_dir);

% Preallocate data array
sweep_results = zeros(length(var_vec), total_runs, 3);




%% Run solvers and save intermediate results

% loop over sweep values
disp('=== Begin solving ===')
for i = 1:length(var_vec)

    % load input sweep parameter
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

        % clear memory
        clear V omega G_stats;

    end

end
disp('=== End solving ===')




%% Save results

% Save data to dir
save([out_dir 'results.mat']);

% Plot results 
TG_plot_sweep_1d(sweep_results, var_str, var_vec, out_dir)




disp('Done!')

end





