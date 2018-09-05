function Params = TG_main_parameters
%myFun - Description
%
% Syntax: Params = myFun(input)
%
% Long description

echo on;

Params.N = 100;  % Total number of neurons
Params.f = 0.5;  % Proportion of excitatory nodes 
Params.tau = 1/4;  % Decay rate of system

Params.sigmae = 5;  % excitatory standard deviation
Params.sigmai = 1;  % inhibitory standard deviation

Params.mue = 0;  % excitatory mean
Params.mui = -Params.mue;  % inhibitory mean

Params.frac_EV_TG = 0.5;  % fraction of eigenvalues used for TG

Params.t_min = 0;
Params.t_max = 10;
Params.t_step = 0.1;

Params.seed = 0;

echo off;

end