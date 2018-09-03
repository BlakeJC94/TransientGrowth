function TG_plot_eigvals(N, f, tau, sigma_vec, mu_vec, omega)
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
sigmae = sigma_vec(1);  % excitatory standard deviation
sigmai = sigma_vec(2);  % inhibitory standard deviation

mue = mu_vec(1);  % excitatory mean
mui = mu_vec(2);  % inhibitory mean



%% Calulate circle radius and plot vecs
sigmaeff = sqrt(f*sigmae^2 + (1-f)*sigmai^2);
theta = linspace(0,2*pi,360);



%% Plot options 
figopt.xlabelstr = '$\omega_r$';
figopt.ylabelstr = '$\omega_i$';
figopt.titlestr = ['$N$=' num2str(N) ',~$f$=' num2str(f) ... 
    '$,~\sigma_{e}$=' num2str(sigmae),'$,~\sigma_{i}$=',num2str(sigmai)];

figopt.xmin = -2*sigmaeff - 1/tau + 1;
figopt.xmax = 1;
figopt.ymin = -1*sigmaeff - 1/(2*tau);
figopt.ymax =  1*sigmaeff + 1/(2*tau);
figopt.fitfactor = 1.3;



%% Plot figure 
figure; 
hold on;
% plot circle enclosing eigenspectrum in solid black
plot(-1/tau + sigmaeff*cos(theta), sigmaeff*sin(theta),'-k','linewidth',1)
% plot dashed grey line from (0, 10) to (0 ,-10) 
plot([0 0],[-10 10],'--','color',[0.5 0.5 0.5])
% plot eigenvalues of J
plot(omega,'o','markersize',4); 
hold off;

box on;

xlabel(figopt.xlabelstr,'interpreter','latex');
ylabel(figopt.ylabelstr,'interpreter','latex');
title(figopt.titlestr,'interpreter','latex');

axis([figopt.xmin figopt.xmax figopt.ymin figopt.ymax]*figopt.fitfactor)
axis square



end


