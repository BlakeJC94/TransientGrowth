function Params = TG_parameters
%TG_main_parameters - Stores and outputs default parameters for TG_main
%
% Syntax: Params = TG_main_parameters
%
% Values specified here will be used as the default set of parameters during 
%   demo mode (TG_main or TG_main(0)) and the default set for 1d sweeps
%
% run display(Params) after Params = TG_parameters to print to CLI




N = 400;  % Total number of neurons
f = 0.5;  % Proportion of excitatory nodes 
tau = 1/40;  % Decay rate of system

sigmae = 1/(sqrt(N)*tau);  % excitatory standard deviation
alpha = 1/5;
sigmai = alpha*sigmae;  % inhibitory standard deviation

mue = 0;  % excitatory mean
mui = -mue;  % inhibitory mean

frac_EV_TG = 0.5;  % fraction of eigenvalues used for TG

R = sqrt(N*(f*sigmae^2 + (1-f)*sigmai^2));

t_min = 0;
t_max = 2/(abs(1/tau - R));
t_step = 100;

seed = 'shuffle';




% wrap paramters
Params.N = N;  % Total number of neurons
Params.f = f;  % Proportion of excitatory nodes 
Params.tau = tau;  % Decay rate of system

Params.sigmae = sigmae;  % excitatory standard deviation
Params.sigmai = sigmai;  % inhibitory standard deviation

Params.mue = mue;  % excitatory mean
Params.mui = mui;  % inhibitory mean

Params.frac_EV_TG = frac_EV_TG;  % fraction of eigenvalues used for TG

Params.t_min = t_min;
Params.t_max = t_max;
Params.t_step = t_step;

Params.seed = seed;


end