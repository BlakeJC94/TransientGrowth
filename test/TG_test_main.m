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
plot_eigvals = 0;
plot_G_vec = 1;
% rng('shuffle');  % randomises results
%=========================================================================
% Variable declaration
N = 1000;  % Total number of neurons
f = 0.5;  % Proportion of excitatory nodes 
tau = 1/4;  % Decay rate of system

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
    TG_get_eig_matrix(N, f, tau, sigma_vec, mu_vec, frac_EV_TG, seed);


if plot_eigvals == 1
    % calculate radius of circular eigenspectrum
    sigmaeff = sqrt(N*(f*sigmaev^2/N + (1-f)*sigmaiv^2/N));
    
    %---------------------plot Eigenvalues--------------------
    % plot options
    fig1opt.xlabelstr = '$\omega_r$';
    fig1opt.ylabelstr = '$\omega_i$';
    fig1opt.titlestr = ['$N$=' num2str(N) ',~$f$=' num2str(f) ... 
        '$,~\sigma_{e}$=' num2str(sigmaev),'$,~\sigma_{i}$=',num2str(sigmaiv)];
    
    fig1opt.xmin = -2*sigmaeff+1;
    fig1opt.xmax = 1;
    fig1opt.ymin = -sigmaeff;
    fig1opt.ymax = sigmaeff;
    fig1opt.fitfactor = 1.3;
    
    % plot figure
    figure; hold on;
    
    theta = linspace(0,2*pi,360);
    
    % plot circle enclosing eigenspectrum in solid black
    plot(-1/tau + sigmaeff*cos(theta), sigmaeff*sin(theta),'-k','linewidth',1)
    % plot dashed grey line from (0, 10) to (0 ,-10) 
    plot([0 0],[-10 10],'--','color',[0.5 0.5 0.5])
    % plot eigenvalues of J
    plot(omega,'o','markersize',4); 
    
    box on;
    
    xlabel(fig1opt.xlabelstr,'interpreter','latex');
    ylabel(fig1opt.ylabelstr,'interpreter','latex');
    title(fig1opt.titlestr,'interpreter','latex');
    
    %axis([fig1opt.xmin fig1opt.xmax fig1opt.ymin fig1opt.ymax]*fig1opt.fitfactor)
    axis square
end



t_vec = t_min:t_step:t_max;
G_vec = TG_get_max_growth(t_vec, V, omega, frac_EV_TG);


[G_max, G_max_ind] = max(G_vec);
t_opt = t_vec(G_max_ind);
G_init_slope = (G_vec(2) - G_vec(1))/(t_vec(2) - t_vec(1));



if plot_G_vec == 1
    % plot options 
    fig2opt.xlabelstr = '$t$';
    fig2opt.ylabelstr = ['$G(t):=\displaystyle\max_{{\bf u(t=0)}}' ...
        '\frac{{\bf u}(t)}{{\bf u}(t=0)}$'];
    
    % plot figure
    figure;
    plot(t_vec, G_vec, '--')
    xlabel(fig2opt.xlabelstr,'interpreter','latex');
    ylabel(fig2opt.ylabelstr,'interpreter','latex');
    
    
    hold on;
    plot([t_opt, t_opt], [G_max, 0], 'k--');
    text(t_opt+0.1, 0.1, ['$t_{opt} = $' num2str(t_opt)], 'Interpreter', 'latex');
    text(t_opt+0.1, G_max+0.1, ['$G_{max} = $' num2str(G_max)], 'Interpreter', 'latex');
    text(t_vec(1)+0.1, G_vec(1), ['$\frac{dG}{dt}|_{t=0} = $' num2str(G_init_slope)], 'Interpreter', 'latex');
    hold off;
end

