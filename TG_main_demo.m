function TG_main_demo(var_str_cell, var_vec)
%TG_main - Main function used to run demos (1 run)
%
% Syntax: TG_main_demo(var_str_cell, var_vec)
% 
% To run for default values, run TG_main_demo with no input arguements. To 
%   change a selection of parameters, set var_str_cell = {'var1', 'var2'} and 
%   var_val = [val1, val2]. Length of inputs must be equal.
% 
% e.g. >> TG_main(0, {'N', 'tau'}, [150, 1/5]) 
%
% Input: 
%   var_str_cell : cell of strings of variable names to change
%   var_vec : vector of values
% 
% Data saved to 'output/TG_main_demo/'


%% Initialise

addpath(genpath('src'))
close all;

% check input
if nargin == 1
    error('Requested parameter change, but no value given!')
elseif nargin == 0
    var_str_cell = [];
    var_vec = [];
end
if length(var_str_cell) ~= length(var_vec)
    error('Length of arguments must be equal')
end

% Set output directory 
out_dir = 'output/TG_main_demo/';




%% Load parameters 

% load default parameters
Params = TG_parameters;
display(Params);

% change parameters if requested
if ~isempty(var_str_cell)

    % Throw error if multiple values given for var_str_cell
    if length(var_str_cell) ~= length(var_vec)
        error('Parameter name value mismatch in demo mode!')
    end

    disp('Parameters modified')
    for i = 1:length(var_str_cell)
        eval(['Params.' var_str_cell{i} ' = ' num2str(var_vec(i)) ';']);
    end
end

% Save input to file
TG_write_input(Params, out_dir);




%% Solve system 

% Get eigenvalues of J
[V, omega, J_mat] = TG_get_eig_matrix(Params);  %#ok<ASGLU>

% Get max growth
[G_vec, G_stats] = TG_get_max_growth(Params, V, omega);




%% Save results

% Write G_stats to file
TG_write_output(G_stats, out_dir);

% Plot eigenvalues
TG_plot_eigvals(Params, omega, out_dir);

% Plot max growth
TG_plot_max_growth(Params, G_vec, G_stats, out_dir);

% Save data to dir
save([out_dir 'results.mat']);




disp('Done!')

end





