% Transient growth calcualtions for Neural Networks: Matrices of the form
% given by Rajan and Abbott (PRL - 2006)
clear; %clc;
tau = -4; f = 0.5;  
sigmaev = 5; %excitatory variance
sigmaiv = 1; %inhibitory variance
Ne=200; Ni=200; N=Ne+Ni; mui=0; mue=-mui; 
%=========================================================================

seed=1; 
% rng(seed);
rng('shuffle');

sigmai1=sigmaiv/sqrt(N);
sigmae1=sigmaev/sqrt(N);
Ai = mui + sigmai1*randn(N,ceil(N/2));
Ae = mue + sigmae1*randn(N,floor(N/2));
J = [Ae Ai];
Tau = eye(size(J))*tau;

[V,omega] = eig(J+Tau);

omega = diag(omega);
sigmaeff = sqrt(N*(f*sigmae1^2+(1-f)*sigmai1^2))  ; % radius of circular eigenspectrum

%---------------------plot Eigenvalues--------------------
figure; hold on;
theta = linspace(0,2*pi,360);
plot(tau + sigmaeff*cos(theta), sigmaeff*sin(theta),'-k','linewidth',1)
plot([0 0],[-10 10],'--','color',[0.5 0.5 0.5])
plot(omega,'o','markersize',4); 
box on;
xlabel('$\omega_r$','interpreter','latex');ylabel('$\omega_i$','interpreter','latex');
title(['$N$=',num2str(N),',~$f$=',num2str(f),'$,~\sigma_{e}$=',num2str(sigmaev),'$,~\sigma_{i}$=',num2str(sigmaiv)], ...
    'interpreter','latex');
axis([-2*sigmaeff+1 1 -sigmaeff sigmaeff]*1.3)
axis square
%%
%-----------------ops05t-- Get rid of spurious+smaller EignValues for TG
frac_EV_TG = 0.5; % fraction of eigenvalues used for TG
n_EG = floor(N*frac_EV_TG);

[~,ii_omega_TG] = sort(real(omega),'descend');
ii_omega_TG = ii_omega_TG(1:n_EG);

%----------------- ops1ch--Get the corresponding EigenFunctions for TG
M = length(ii_omega_TG); % = n_EG = the number of Eigenvalues used for expansion
V_u = zeros(N,M);
for jj = 1:M
    V_u(:,jj) = V(1:N,ii_omega_TG(jj));
end
% clear V 
%---------------------plot EigenFunctions--------------------
% figure;
% plot(abs(V_u(:,1)),'.-')


%-----------------ops31ch--'u': Integrate A components
A_u = zeros(M,M);
for ii_col = 1:M % it's awaste of resources because A_u is symm!, but readability
    for ii_row = 1:M
        A_u(ii_row,ii_col) = sum(conj(V_u(:,ii_row)).*V_u(:,ii_col));  
    end
end
A = A_u;
clear A_u
%==========================================================================
% [[ help eig ==> [V,D] = eig(A) returns A*V = V*D. ]]
% We want SVD(F*LamX*inv(F)), where, A = F'F. To find F --->
% Now, A V =  V lam => A = V lam V'; Note V' = inv(V) --> A is symm postive
% F = V lam V', because F'F=(V lam V')'(V lam V)=V lam' V' V lam V' =
% V lam^2 V' = A.
% Also, LamX = exp(-i omega t) --> See notes
%==========================================================================
[V_a,lam_a] = eig(A);
F_a = V_a*sqrt(lam_a)*V_a';

t = 0:0.1:5;
G = zeros(length(t),1);
for ii = 1:length(t)  %  t(non-dim) varies [1:Re/10] and peaks @ t/Re=0.08[Thesis p.24]
    LamX = diag(exp(omega(ii_omega_TG)*t(ii)));
%     LamX = diag(exp(-i*omega(ii_omega_TG)*300));
    [U_opt,S,V_opt] = svd(F_a*LamX*inv(F_a)); %[U,S,V] = svd(X) produces X*V = U*S.
    G(ii) = S(1,1)^2;
end
figure;
plot(t, G, '--')
xlabel('$t$','interpreter','latex');
ylabel('$G(t):=\displaystyle\max_{{\bf u(t=0)}}\frac{{\bf u}(t)}{{\bf u}(t=0)}$','interpreter','latex');






