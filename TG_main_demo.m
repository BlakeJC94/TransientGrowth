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
[V, omega, J_mat] = TG_get_eig_matrix(Params); 

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



%% Solve ODE system

% specify initial condition
u_init = 0.0001*randn(Params.N, 1);

t_vec = linspace(Params.t_min, Params.t_max, Params.t_step);
t_dif = abs(t_vec(2) - t_vec(1));

% preallocate memory (index by (space, time))
u_mat = zeros(Params.N, length(t_vec));

% use forward euler solve
u_vec = u_init;
u_mat(:, 1) = u_vec;
for ind = 2:length(t_vec)

    % Solve next step
    u_vec = (eye(Params.N) + t_dif*J_mat) * u_vec;
    
    % save result to u_mat
    u_mat(:, ind) = u_vec;
    
end

figure;
imagesc(1:Params.N, t_vec, abs(u_mat')); 
set(gca,'YDir','normal');
xlabel('Spatial index');
ylabel('Time');
% text(1, G_stats.t_opt, '$$t_{opt}$$', 'interpreter', 'latex');
colormap('jet');
colorbar;

hold on;
plot(1:Params.N, (G_stats.t_opt)*ones(1,Params.N), 'k--')
hold off;

plot_export_fig(0, [out_dir 'figures/plot_spacetime_heat'], 14, 7/5, 18);


disp('Done!')

end





