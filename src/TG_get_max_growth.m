function G_vec = TG_get_max_growth(t_vec, V, omega, frac_EV_TG)
%TG_get_max_growth - Calculates maximum growth given time vector and 
%   eigenvals/vecs from J and fraction of eigenvalues to keep
%
% Syntax: G_vec = TG_get_max_growth(t_vec, V, omega, frac_EV_TG)
%
% Long description



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


