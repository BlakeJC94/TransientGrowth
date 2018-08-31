% Transient growth calcualtions for Neural Networks: Matrices of the form
% given by Rajan and Abbott (PRL - 2006)
% modified by Blake (Aug 2018)
%=========================================================================
% Transient growth in system given by 
%   \frac{du_i}{dt} = -\frac{u_i}{tau} + \Sum_j W_{ij} \phi(x_j)
%=========================================================================
clear; %clc;
close all;
seed=1; rng(seed);  % makes code repeatable
% rng('shuffle');  % randomises results
%=========================================================================
% Variable declaration
tau = 1/4;  % Decay rate of syste,
f = 0.5;  % Proportion of excitatory nodes 

sigmaev = 5;  % excitatory standard deviation
sigmaiv = 1;  % inhibitory standard deviation

mue = 0;  % excitatory mean
mui = -mue;  % inhibitory mean

N = 400;  % Total number of neurons
Ne = round(f*N); Ni = N - Ne; % Number of excitatory and inhibitory neurons
%=========================================================================


% Create random matrix J = W - I/tau such that with first Ne columns of W are 
% distributed as N(mue, sigmaev^2/N) and the next of Ni columns of W  
% distributed as N(mui, sigmai^2/N)
sigmae1=sigmaev/sqrt(N);
sigmai1=sigmaiv/sqrt(N);

We = mue + sigmae1*randn(N,Ne);
Wi = mui + sigmai1*randn(N,Ni);

W = [We Wi];
Tau = -eye(size(W))/tau;

J = W + Tau;

% obtain eigenvalues and eigenvectors of J
[V,omega] = eig(J);
omega = diag(omega);

% calculate radius of circular eigenspectrum
sigmaeff = sqrt(N*(f*sigmae1^2+(1-f)*sigmai1^2));

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

axis([fig1opt.xmin fig1opt.xmax fig1opt.ymin fig1opt.ymax]*fig1opt.fitfactor)
axis square
%%sigmaeff
%-----------------ops05t-- Get rid of spurious+smaller EignValues for TG
frac_EV_TG = 0.5; % fraction of eigenvalues used for TG
n_EG = floor(N*frac_EV_TG);

% Store indicies of the n_EG eigenvalues with largest real part
[~,ii_omega_TG] = sort(real(omega),'descend');
ii_omega_TG = ii_omega_TG(1:n_EG);

%----------------- ops1ch--Get the corresponding EigenFunctions for TG
M = length(ii_omega_TG); % = n_EG = the number of Eigenvalues used for expansion
% store eigenvectors corresponding to selected eigenvalues above in V_u
V_u = zeros(N,M);
for jj = 1:M
    V_u(:,jj) = V(1:N,ii_omega_TG(jj));
end
% clear V 
%---------------------plot EigenFunctions--------------------
% figure;
% plot(abs(V_u(:,1)),'.-')


%-----------------ops31ch--'u': Integrate A components
% see notes: construct A matrix such that 
%   A_{i,j} = \langle u_i, u_j \rangle = \bar(u_i) \cdot u_j,
% where u_i is the i'th eigenvector (indexing is different here, clearly)
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
[V_a, lam_a] = eig(A);
F_a = V_a*sqrt(lam_a)*V_a';
%disp(diag(V_a * sqrt(lam_a) * (V_a)' *)
F_a_inv =  V_a * diag(1./diag(sqrt(lam_a))) * V_a';

t = 0:0.1:5;
G = zeros(length(t),1);
for ii = 1:length(t)  %  t(non-dim) varies [1:Re/10] and peaks @ t/Re=0.08[Thesis p.24]
    LamX = diag(exp(omega(ii_omega_TG)*t(ii)));

    [U_opt,S,V_opt] = svd(F_a*LamX*F_a_inv); %[U,S,V] = svd(X) produces X*V = U*S.
    G(ii) = S(1,1)^2;
end

% plot options 
fig2opt.xlabelstr = '$t$';
fig2opt.ylabelstr = ['$G(t):=\displaystyle\max_{{\bf u(t=0)}}' ...
    '\frac{{\bf u}(t)}{{\bf u}(t=0)}$'];

% plot figure
figure;
plot(t, G, '--')
xlabel(fig2opt.xlabelstr,'interpreter','latex');
ylabel(fig2opt.ylabelstr,'interpreter','latex');


%%%%%%%%
% TODO: 
%   - Extract statistics from G
%   - Write basic parameter sweeps for G over tau
%   - Sweep over N
%   - 2D sweep over alpha and f
%   - Explore methods of looking at the initial slope
%     - Analytically, paper from Schmidt is a good place to start (p137)


[G_max, G_max_ind] = max(G);
t_max = t(G_max_ind);

hold on;
plot([t_max, t_max], [G_max, 0], 'k--');
text(t_max+0.1, 0.1, ['$t_{max} = $' num2str(t_max)], 'Interpreter', 'latex');
text(t_max+0.1, G_max+0.1, ['$G_{max} = $' num2str(G_max)], 'Interpreter', 'latex');
hold off;
