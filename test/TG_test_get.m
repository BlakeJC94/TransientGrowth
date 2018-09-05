function TG_test_get
% TG_test_get - Tests solvers

addpath(genpath('src'))
close all;

%% load verified data
load('test/data/TG_test_get/results.mat')

% Varible declaration
Params.N = 400;  % Total number of neurons
Params.f = 0.5;  % Proportion of excitatory nodes 
Params.tau = 1/4;  % Decay rate of system

Params.sigmae = 5;  % excitatory standard deviation
Params.sigmai = 1;  % inhibitory standard deviation

Params.mue = 0;  % excitatory mean
Params.mui = -Params.mue;  % inhibitory mean

Params.frac_EV_TG = 0.5;  % fraction of eigenvalues used for TG

Params.t_min = 0;
Params.t_max = 6;
Params.t_step = 0.1;

Params.seed = 0;



%% Test TG_get_eig_matrix()
% Generate J matrix
[V_test, omega_test, J_mat_test] = TG_get_eig_matrix(Params);

% Compare output of TG_get_eig_matrix() to verified results
V_check_flag = isequal(V_test, V);
omega_check_flag = isequal(omega_test, omega);
J_mat_check_flag = isequal(J_mat_test, J_mat);

TG_get_eig_matrix_result = ...
    min([V_check_flag, omega_check_flag, J_mat_check_flag]);

if TG_get_eig_matrix_result
    disp('TG_test_get, TG_get_eig_matrix : Pass!');
else
    warning('TG_test_get, TG_get_eig_matrix : Fail');
end



%% Test TG_get_max_growth()
% Generate G_vec
[G_vec_test, G_stats_test] = TG_get_max_growth(Params, V, omega);

% Compare output of TG_get_eig_matrix() to verified results
G_vec_check_flag = isequal(G_vec_test, G_vec);
G_stats_check_flag = isequal(G_stats_test, G_stats);

TG_get_max_growth_result = ...
    min([G_vec_check_flag, G_stats_check_flag]);

if TG_get_max_growth_result
    disp('TG_test_get, TG_get_max_growth : Pass!');
else
    warning('TG_test_get, TG_get_max_growth : Fail');
end



%% Test complete
TG_get_result = min([TG_get_eig_matrix_result, TG_get_max_growth_result]);
if TG_get_result
    disp('TG_test_get : Pass!');
else 
    warning('TG_test_get : Fail');
end



end