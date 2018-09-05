% Variable declaration
N = 100;  % Total number of neurons
f = 0.5;  % Proportion of excitatory nodes 
tau = 1/4;  % Decay rate of system

sigmae = 5;  % excitatory standard deviation
sigmai = 1;  % inhibitory standard deviation

mue = 0;  % excitatory mean
mui = -mue;  % inhibitory mean

frac_EV_TG = 0.5;  % fraction of eigenvalues used for TG

t_min = 0;
t_max = 10;
t_step = 0.1;

seed = 0;
