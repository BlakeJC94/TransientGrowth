function TG_plot_eigvals(Params, omega, out_dir)
%TG_plot_eigvals - Plots eigenvalues and radius of density in argand plane
%
% Syntax: TG_plot_eigvals(N, f, tau, sigma_vec, mu_vec)
%
% Simple plotting component for looking at eigenvalues of J overlayed on the 
% circle enclosing the spectral density, returns figure
% 
% Inputs:
%   N = Total number of neurons
%   tau = Decay rate of syste,
%   f = Proportion of excitatory nodes 
%   sigma_vec(1), signma_vec(2) = excitatory/inhibitory standard deviation
%   mu_vec(1), mu_vec(2) = excitatory/inhibitory mean
%   omega = eigenvalues of J (from TG_get_eig_matrix())
    


%% Load parameters
N = Params.N;  % Total number of neurons
f = Params.f;  % Proportion of excitatory nodes 
tau = Params.tau;  % Decay rate of system

sigmae = Params.sigmae;  % excitatory standard deviation
sigmai = Params.sigmai;  % inhibitory standard deviation

% mue = Params.mue;  % excitatory mean
% mui = Params.mui;  % inhibitory mean



%% Calulate circle radius and plot vecs
sigmaeff = sqrt(N*(f*sigmae^2 + (1-f)*sigmai^2));
theta = linspace(0,2*pi,360);



%% Plot options 
figopt.xlabelstr = '$\omega_r$';
figopt.ylabelstr = '$\omega_i$';
figopt.titlestr = ['$N$=' num2str(N) ',~$f$=' num2str(f) ... 
    '$,~\sigma_{e}$=' num2str(sigmae),'$,~\sigma_{i}$=',num2str(sigmai)];

figopt.xmin = -2*sigmaeff - 1/tau + 1 + max(0, max(real(omega)));
figopt.xmax = 1 + max(0, max(real(omega)));
figopt.ymin = -1*sigmaeff - 1/(2*tau);
figopt.ymax =  1*sigmaeff + 1/(2*tau);
figopt.fitfactor = 1.3;



%% Plot figure 
figure; 
hold on;
% plot circle enclosing eigenspectrum in solid black
plot(-1/tau + sigmaeff*cos(theta), sigmaeff*sin(theta),'-k','linewidth',1)
% plot dashed grey line from (0, ymin) to (0, ymax) 
plot([0, 0],[figopt.ymin figopt.ymax],'--', 'color', [0.5 0.5 0.5])
% plot dashed grey line from (xmin, 0) to (xmax, 0) 
plot([figopt.xmin figopt.xmax], [0,0],'--', 'color', [0.5 0.5 0.5])
% plot eigenvalues of J
plot(omega,'ro','markersize',4); 
hold off;

box on;

xlabel(figopt.xlabelstr,'interpreter','latex');
ylabel(figopt.ylabelstr,'interpreter','latex');
title(figopt.titlestr,'interpreter','latex');

axis([figopt.xmin figopt.xmax figopt.ymin figopt.ymax]*figopt.fitfactor)
axis square

plot_export_fig(-1, [out_dir 'figures/plot_eigvals'], 14, 7/5, 18);



end


