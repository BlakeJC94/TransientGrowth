function [G_vec, G_stats] = TG_get_max_growth(Params, V, omega)
%TG_get_max_growth - Calculates maximum growth and key stats given time vector 
%   and eigenvals/vecs from J and fraction of eigenvalues to keep 
%
% Syntax: G_vec = TG_get_max_growth(t_vec, V, omega, frac_EV_TG)
%
% Transient growth in system given by 
%     \frac{du_i}{dt} = -\frac{u_i}{tau} + \Sum_j W_{i,j} u_j %TODO: UPDATE THIS AND FOLLOWING CODE
% such that W_{i,j} ~ N(\mu_e, \sigma_e^2/N) for j =< round(f*N) 
%   and W_{i,j} ~ N(\mu_i, \sigma_i^2/N) for j > round(f*N) 
% 
% 1. Extract fraction of eigenstuff with leargest real part
% 2. Construct F Lambda F_inv given eigenstuff from J
% 3. Calculate operator norm of F Lambda F_inv
% 
% Input:
%   t_vec = vector of values of t
%   V = eigenvectors of J (from TG_get_eig_matrix)
%   omega = eigenvalues of J (from TG_get_eig_matrix)
%   frac_EV_TG = fraction of eigenvalues to use
% 
% Output: 
%   G_vec = vector of values of G for intput t
%   G_stats = struct containing key stats for G
%     G_max = maximum amplification, (local max found using max)
%     t_opt = time at which max occurs
%     G_init_slope = slope of G at t=min(t)^+ (found using 1st order diff)
% 
% 

frac_EV_TG = Params.frac_EV_TG;
t_vec = linspace(Params.t_min, Params.t_max, Params.t_step);

%% Eliminate spurious/smaller eigenvalues
N = length(omega);
n_EG = floor(N*frac_EV_TG);
% store indices of the n_EG eigenvalues with largest real part
% (most likely to casue instabilities)
[~, ii_omega_TG] = sort(real(omega), 'descend');
ii_omega_TG = ii_omega_TG(1:n_EG);
% get corresponding eigenvectors 
M = length(ii_omega_TG);
V_u = zeros(N, M);
for index = 1:M
    V_u(:, index) = V(:, ii_omega_TG(index));
end



%% Construct A 
% A_{i,j} = \langle v_i, v_j \rangle = \bar(v_i) \cdot v_j,
% where v_i is the i'th eigenvector (indexing is different here, clearly)
A_u = zeros(M, M);
% Calculate first diag entry
A_u(1, 1) = sum(conj(V_u(:, 1)).*V_u(:, 1));
for row = 2:M 
    % Calculate lower off-diag entries
    for col = 1:row-1
        A_u(row, col) = sum(conj(V_u(:, row)).*V_u(:, col));
    end
    % Calculate diag entries
    A_u(row, row) = sum(conj(V_u(:,row)).*V_u(:,row));
end
% fill in upper off-diag entries
A = A_u + tril(A_u,-1)';



%% Construct F such that A = F'* F 
% ORIGINAL METHOD
% [V_a, lam_a] = eig(A);
% F_a = V_a*sqrt(lam_a)*V_a';
% F_a_inv =  V_a * diag(1./diag(sqrt(lam_a))) * V_a';
% NEW METHOD, slightly faster
F_a = chol(A);
F_a_inv = F_a\eye(size(F_a));



%% Construct Lambda and G_vec
G_vec = zeros(length(t_vec),1);
for index = 1:length(t_vec)
    Lambda = diag(exp(omega(ii_omega_TG)*t_vec(index)));
    % Matrix 2-norm equivilent to largest singular value
    G_vec(index) = norm(F_a*Lambda*F_a_inv)^2;
end



%% Calculate key statistics of G_vec
[G_max, G_max_ind] = max(G_vec); 
t_opt = t_vec(G_max_ind);
G_init_slope = (G_vec(2) - G_vec(1))/(t_vec(2) - t_vec(1));
% check if peak is truly a local max
if abs(G_max_ind - length(G_vec)) < 10
    warning('Peak occurs in last 10 steps, increase t_max')
end
% export stats
G_stats.G_max = G_max;
G_stats.t_opt = t_opt;
G_stats.G_init_slope = G_init_slope;















% %-----------------ops05t-- Get rid of spurious+smaller EignValues for TG
% n_EG = floor(size(omega,1)*frac_EV_TG);

% % Store indicies of the n_EG eigenvalues with largest real part
% [~,ii_omega_TG] = sort(real(omega),'descend');
% ii_omega_TG = ii_omega_TG(1:n_EG);



% [V_a, lam_a] = eig(A_mat);
% F_a = V_a*sqrt(lam_a)*V_a';
% %disp(diag(V_a * sqrt(lam_a) * (V_a)' *)
% F_a_inv =  V_a * diag(1./diag(sqrt(lam_a))) * V_a';

% G_vec = zeros(length(t_vec),1);
% for ii = 1:length(t_vec)  %  t(non-dim) varies [1:Re/10] and peaks @ t/Re=0.08[Thesis p.24]
%     LamX = diag(exp(omega(ii_omega_TG)*t_vec(ii)));

%     [~, S, ~] = svd(F_a*LamX*F_a_inv); %[U,S,V] = svd(X) produces X*V = U*S.
%     G_vec(ii) = S(1,1)^2;
% end

% [~, G_max_ind] = max(G_vec);

% if abs(G_max_ind - length(t_vec)) < 1
%     warning('G_max detected to be close to t_max, check plots and increase t_max')
% end

% end


