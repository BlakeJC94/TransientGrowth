function TG_test_get_results
% TG_test_get_results - Generates data for test, MATLAB R2018a

addpath(genpath('src'));

% Varible declaration
N = 400;  % Total number of neurons
f = 0.5;  % Proportion of excitatory nodes 
tau = 1/4;  % Decay rate of system

sigmae = 5;  % excitatory standard deviation
sigmai = 1;  % inhibitory standard deviation

mue = 0;  % excitatory mean
mui = -mue;  % inhibitory mean

frac_EV_TG = 0.5;  % fraction of eigenvalues used for TG

t_min = 0;
t_max = 6;
t_step = 0.1;

seed = 0;

% set up parameter inputs
sigma_vec = [sigmae, sigmai];
mu_vec = [mue, mui];
t_vec = t_min:t_step:t_max;



% Generate J matrix
[V, omega, J_mat] = ...
    TG_get_eig_matrix(N, f, tau, sigma_vec, mu_vec, seed);

% Generate G_vec
[G_vec, G_stats] = TG_get_max_growth(t_vec, V, omega, frac_EV_TG);



%% save verified data
clear N f tau sigmae sigmai sigma_vec mue mui mu_vec frac_EV_TG t_min t_max ...
    t_step t_vec seed;
save('test/data/TG_test_get/results.mat')


    
end