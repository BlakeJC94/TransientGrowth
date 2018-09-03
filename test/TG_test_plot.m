function TG_test_plot
%TG_test_plot - Tests plotting functions 
% Verified results located in 'output/test/TG_test_plot/verified/'

addpath(genpath('src'))
close all;


%% Load data
% Using original set of paramter values:
% Varible declaration
% N = 400;  % Total number of neurons
% f = 0.5;  % Proportion of excitatory nodes 
% tau = 1/4;  % Decay rate of system
% sigmae = 5;  % excitatory standard deviation
% sigmai = 1;  % inhibitory standard deviation
% mue = 0;  % excitatory mean
% mui = -mue;  % inhibitory mean
% frac_EV_TG = 0.5;  % fraction of eigenvalues used for TG
% t_min = 0;
% t_max = 6;
% t_step = 0.1;
% seed = 0;
load('test/data/TG_test_plot/results.mat')
fig_dir = 'output/test/TG_test_plot/';

%% Test TG_plot_eigvals(N, f, tau, sigma_vec, mu_vec, omega)
TG_plot_eigvals(N, f, tau, sigma_vec, mu_vec, omega);
SaveAsPngEpsAndFig(-1, [fig_dir 'TG_plot_eigvals_results'], 14, 7/5, 18);

%% Test TG_plot_max_growth(G_vec, t_vec, G_stats)
TG_plot_max_growth(G_vec, t_vec, G_stats)
SaveAsPngEpsAndFig(-1, [fig_dir 'TG_plot_max_growth_results'], 14, 7/5, 18);

%% Test complete
disp('TG_test_plot : Pass!')

end




% function TG_test_plot_gen_data
% % Run if data file is missing
% addpath(genpath('src'))

% % Varible declaration
% N = 400;  % Total number of neurons
% f = 0.5;  % Proportion of excitatory nodes 
% tau = 1/4;  % Decay rate of system

% sigmae = 5;  % excitatory standard deviation
% sigmai = 1;  % inhibitory standard deviation

% mue = 0;  % excitatory mean
% mui = -mue;  % inhibitory mean

% frac_EV_TG = 0.5;  % fraction of eigenvalues used for TG

% t_min = 0;
% t_max = 6;
% t_step = 0.1;

% seed = 0;


% % set up parameter inputs
% sigma_vec = [sigmae, sigmai];
% mu_vec = [mue, mui];
% t_vec = t_min:t_step:t_max;


% % Generate J matrix
% [V, omega, J_mat] = ...
%     TG_get_eig_matrix(N, f, tau, sigma_vec, mu_vec, seed);

% % Generate G_vec
% [G_vec, G_stats] = TG_get_max_growth(t_vec, V, omega, frac_EV_TG);

% % Save data to file
% save('test/data/TG_test_plot/results.mat')


% end