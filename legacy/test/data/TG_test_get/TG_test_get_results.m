function TG_test_get_results
% TG_test_get_results - Generates data for test, MATLAB R2018a

addpath(genpath('src'));

% Varible declaration
Params.N = 400;  % Total number of neurons
Params.f = 0.5;  % Proportion of excitatory nodes 
Params.tau = 1/4;  % Decay rate of system

Params.sigmae = 5;  % excitatory standard deviation
Params.sigmai = 1;  % inhibitory standard deviation

Params.mue = 0;  % excitatory mean
Params.mui = -mue;  % inhibitory mean

Params.frac_EV_TG = 0.5;  % fraction of eigenvalues used for TG

Params.t_min = 0;
Params.t_max = 6;
Params.t_step = 0.1;

Params.seed = 0;



% Generate J matrix
[V, omega, J_mat] = TG_get_eig_matrix(Params);

% Generate G_vec
[G_vec, G_stats] = TG_get_max_growth(Params, V, omega);



%% save verified data
save('test/data/TG_test_get/results.mat')


    
end