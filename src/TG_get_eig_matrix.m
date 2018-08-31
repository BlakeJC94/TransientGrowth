function [V, omega, J_mat] = TG_get_eig_matrix(N, f, tau, sigma_vec, mu_vec, frac_EV_TG, seed)
%TG_get_eig_matrix - Constructs the matrix operator J and returns eignevals/vecs
%
% Syntax: [V, omega, J_mat] = 
%     TG_get_eig_matrix(N, f, tau, sigma_vec, mu_vec, frac_EV_TG, seed)
%
% Transient growth in system given by 
%     \frac{du_i}{dt} = -\frac{u_i}{tau} + \Sum_j W_{i,j} u_j %TODO: UPDATE THIS AND FOLLOWING CODE
% such that W_{i,j} ~ N(\mu_e, \sigma_e^2/N) for j =< round(f*N) 
%   and W_{i,j} ~ N(\mu_i, \sigma_i^2/N) for j > round(f*N) 
% 
% 1. Construct J given parameters
% 2. Returns eigenvals/vecs for J
% 
% Input:
%   N = Total number of neurons
%   tau = Decay rate of syste,
%   f = Proportion of excitatory nodes 
%   
%   sigma_vec(1), signma_vec(2) = excitatory/inhibitory standard deviation
%   mu_vec(1), mu_vec(2) = excitatory/inhibitory mean
%   
%   Ne, Ni = Number of excitatory and inhibitory neurons



%% Load parameters 
rng(seed) % put seed = 'shuffle' if random results are desired

sigmae = sigma_vec(1);  % excitatory standard deviation
sigmai = sigma_vec(2);  % inhibitory standard deviation

mue = mu_vec(1);  % excitatory mean
mui = mu_vec(2);  % inhibitory mean

Ne = round(f*N); Ni = N - Ne; % Number of excitatory and inhibitory neurons



%% Create random matrx J = W - I/tau such that 
% first Ne columns of W are distributed as N(mue, sigmaev^2/N) 
% and the next of Ni columns of W distributed as N(mui, sigmai^2/N)
We = mue + (sigmae/sqrt(N))*randn(N, Ne);
Wi = mui + (sigmai/sqrt(N))*randn(N, Ni);

W = [We Wi];
Tau = -eye(size(W))/tau;

J_mat = W + Tau;



%% obtain eigenvalues and eigenvectors of J
[V,omega] = eig(J_mat);
omega = diag(omega);


end


% %-----------------ops05t-- Get rid of spurious+smaller EignValues for TG
% n_EG = floor(N*frac_EV_TG);

% % Store indicies of the n_EG eigenvalues with largest real part
% [~,ii_omega_TG] = sort(real(omega),'descend');
% ii_omega_TG = ii_omega_TG(1:n_EG);

% %----------------- ops1ch--Get the corresponding EigenFunctions for TG
% M = length(ii_omega_TG); % = n_EG = the number of Eigenvalues used for expansion
% % store eigenvectors corresponding to selected eigenvalues above in V_u
% V_u = zeros(N,M);
% for jj = 1:M
%     V_u(:,jj) = V(1:N,ii_omega_TG(jj));
% end

% %-----------------ops31ch--'u': Integrate A components
% % see notes: construct A matrix such that 
% %   A_{i,j} = \langle u_i, u_j \rangle = \bar(u_i) \cdot u_j,
% % where u_i is the i'th eigenvector (indexing is different here, clearly)
% A_u = zeros(M,M);
% A_u(1,1) = sum(conj(V_u(:,1)).*V_u(:,1))/2;
% for ii_row = 2:M
%     for ii_col = 1:ii_row-1
%         A_u(ii_row,ii_col) = sum(conj(V_u(:,ii_row)).*V_u(:,ii_col));
%     end
%     A_u(ii_row,ii_row) = sum(conj(V_u(:,ii_row)).*V_u(:,ii_col))/2;
% end

% A_mat = A_u + A_u';
% clear A_u


% end