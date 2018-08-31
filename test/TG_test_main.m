% Transient growth calcualtions for Neural Networks: Matrices of the form
% given by Rajan and Abbott (PRL - 2006)
% modified by Blake (Aug 2018)
%=========================================================================
% Transient growth in system given by 
%   \frac{du_i}{dt} = -\frac{u_i}{tau} + \Sum_j W_{ij} \phi(x_j)
%=========================================================================
clear; %clc;
close all;
addpath(genpath('../src'))
plot_eigvals = 1;
plot_G_vec = 1;
% rng('shuffle');  % randomises results
%=========================================================================
% Variable declaration
N = 400;  % Total number of neurons
f = 0.5;  % Proportion of excitatory nodes 
tau = 1/8;  % Decay rate of system

sigmaev = 5;  % excitatory standard deviation
sigmaiv = 1;  % inhibitory standard deviation

mue = 0;  % excitatory mean
mui = -mue;  % inhibitory mean

frac_EV_TG = 0.5;  % fraction of eigenvalues used for TG

t_min = 0;
t_max = 10;
t_step = 0.1;

%seed = 'shuffle';
seed = 0;

%=========================================================================

sigma_vec = [sigmaev, sigmaiv];
mu_vec = [mue, mui];


[V, omega, J_mat] = ...
    TG_get_eig_matrix(N, f, tau, sigma_vec, mu_vec, seed);


if plot_eigvals == 1
    TG_plot_eigvals(N, f, tau, sigma_vec, mu_vec, omega);
end



t_vec = t_min:t_step:t_max;
[G_vec, G_stats] = TG_get_max_growth(t_vec, V, omega, frac_EV_TG);


if plot_G_vec == 1
    TG_plot_max_growth(G_vec, t_vec, G_stats)
end

