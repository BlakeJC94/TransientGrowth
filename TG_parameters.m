function Params = TG_parameters
%TG_main_parameters - Stores and outputs default parameters for TG_main
%
% Syntax: Params = TG_main_parameters
%
% Values specified here will be used as the default set of parameters during 
%   demo mode (TG_main or TG_main(0)) and the default set for 1d sweeps
%
% run display(Params) after Params = TG_parameters to print to CLI

Params.N = 400;  % Total number of neurons
Params.f = 0.5;  % Proportion of excitatory nodes 
Params.tau = 1/40;  % Decay rate of system

Params.sigmae = 5;  % excitatory standard deviation
Params.sigmai = 1;  % inhibitory standard deviation

Params.mue = 0;  % excitatory mean
Params.mui = -Params.mue;  % inhibitory mean

Params.frac_EV_TG = 0.5;  % fraction of eigenvalues used for TG

Params.t_min = 0;
R = sqrt(Params.N*(Params.f*Params.sigmae^2 + (1-Params.f)*Params.sigmai^2));
Params.t_max = 3/(abs(1/Params.tau - R));
Params.t_step = 100;

Params.seed = 'shuffle';


end