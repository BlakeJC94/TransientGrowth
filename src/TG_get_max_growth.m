function G_vec = TG_get_max_growth(t_vec, A_mat, omega, frac_EV_TG)
%TG_get_max_growth - Calculates maximum growth given time vec and matrix A
%
% Syntax: G_vec = TG_get_max_growth(t_vec, A_mat)
%
% Long description

%-----------------ops05t-- Get rid of spurious+smaller EignValues for TG
n_EG = floor(size(omega,1)*frac_EV_TG);

% Store indicies of the n_EG eigenvalues with largest real part
[~,ii_omega_TG] = sort(real(omega),'descend');
ii_omega_TG = ii_omega_TG(1:n_EG);



[V_a, lam_a] = eig(A_mat);
F_a = V_a*sqrt(lam_a)*V_a';
%disp(diag(V_a * sqrt(lam_a) * (V_a)' *)
F_a_inv =  V_a * diag(1./diag(sqrt(lam_a))) * V_a';

G_vec = zeros(length(t_vec),1);
for ii = 1:length(t_vec)  %  t(non-dim) varies [1:Re/10] and peaks @ t/Re=0.08[Thesis p.24]
    LamX = diag(exp(omega(ii_omega_TG)*t_vec(ii)));

    [~, S, ~] = svd(F_a*LamX*F_a_inv); %[U,S,V] = svd(X) produces X*V = U*S.
    G_vec(ii) = S(1,1)^2;
end

[~, G_max_ind] = max(G_vec);

if abs(G_max_ind - length(t_vec)) < 1
    warning('G_max detected to be close to t_max, check plots and increase t_max')
end

end


