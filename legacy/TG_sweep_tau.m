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
% rng('shuffle');  % randomises results
%=========================================================================
% Variable declaration
N = 300;  % Total number of neurons
f = 0.5;  % Proportion of excitatory nodes 
%tau = 1/4;  % Decay rate of system

sigmaev = 5;  % excitatory standard deviation
sigmaiv = 1;  % inhibitory standard deviation

mue = 0;  % excitatory mean
mui = -mue;  % inhibitory mean

frac_EV_TG = 0.5;  % fraction of eigenvalues used for TG

t_min = 0;
t_max = 20;
t_step = 0.1;

runs = 20;
tau_vec = 0.125:0.125:1;
%=========================================================================

sigma_vec = [sigmaev, sigmaiv];
mu_vec = [mue, mui];

t_vec = t_min:t_step:t_max;

G_max_mat = zeros(length(tau_vec), runs);
t_opt_mat = zeros(length(tau_vec), runs);
G_init_slope_mat = zeros(length(tau_vec), runs);

for i = 1:length(tau_vec)
    tau = tau_vec(i);
    disp(['  tau = ' num2str(tau)])

    for j = 1:runs
        disp(['    run = ' num2str(j)])

        seed = 1000 * j;


        [A_mat, omega] = TG_get_eig_matrix(N, f, tau, sigma_vec, mu_vec, frac_EV_TG, seed);
        G_vec = TG_get_max_growth(t_vec, A_mat, omega, frac_EV_TG);


        [G_max, G_max_ind] = max(G_vec);
        G_max_mat(i, j) = G_max;
        t_opt_mat(i, j) = t_vec(G_max_ind);
        G_init_slope_mat(i, j) = (G_vec(2) - G_vec(1))/(t_vec(2) - t_vec(1));


        clear A_mat omega G_vec;

    end
    
end

G_max_avg = mean(G_max_mat, 2);
t_opt_avg = mean(t_opt_mat, 2);
G_init_slope_avg = mean(G_init_slope_mat, 2);

ind = sum((t_max - t_opt_mat < 0.5), 2)>0;



fig_dir = '../Figures/tau_sweep/';


figure;
plot(tau_vec, G_max_avg, 'ko-', ...
    tau_vec(ind), G_max_avg(ind), 'r*');
xlabel('N')
title(['G\_max average over ' num2str(runs) ' runs'])
SaveAsPngEpsAndFig(-1, [fig_dir 'g_max_avg'], 14, 7/5, 18);


figure;
plot(tau_vec, t_opt_avg, 'ko-', ...
    tau_vec(ind), t_opt_avg(ind), 'r*');
xlabel('N')
title(['t\_opt average over ' num2str(runs) ' runs'])
SaveAsPngEpsAndFig(-1, [fig_dir 't_opt_avg'], 14, 7/5, 18);


figure;
plot(tau_vec, G_init_slope_avg, 'ko-', ...
    tau_vec(ind), G_init_slope_avg(ind), 'r*');
xlabel('N')
title(['G\_init\_slope average over ' num2str(runs) ' runs'])
SaveAsPngEpsAndFig(-1, [fig_dir 'G_init_slope_avg'], 14, 7/5, 18);
