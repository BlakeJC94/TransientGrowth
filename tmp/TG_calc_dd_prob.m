function prob = TG_calc_dd_prob(Params)
%TG_calc_dd_prob - Calculates probability of diagonal domiance fail numerically
%
% Syntax: output = TG_calc_dd_prob(input)
%
% Calculates probability of the event 
%   {|W_{i,i} - tau^{-1}| < Sum_{j \neq i}|W_{i,j}| for at least 1 i}
%       = {..., for at least 1 i \leq \ceil(N*f)} 
%           \cup {..., for at least 1 i > \ceil(N*f)} 

%% load default parameters
if nargin < 1
    Params = TG_main_parameters;
end
display(Params);


Ne = ceil(Params.f*Params.N); 
Ni = Params.N - Ne;


runs = 100000;

triale_count = 0;
for i = 1:runs
    We = Params.mue + Params.sigmae * randn(1, Ne);
    Wi = Params.mui + Params.sigmai * randn(1, Ni);

    triale = (abs(We(1) - Params.tau^(-1)) ...
        < sum(abs(We(2:end))) + sum(abs(Wi))  );
    
    if triale
        triale_count = triale_count + triale;
    end

end

triali_count = 0;
for i = 1:runs
    We = Params.mue + Params.sigmae * randn(1, Ne);
    Wi = Params.mui + Params.sigmai * randn(1, Ni);

    triali = (abs(Wi(1) - Params.tau^(-1)) ...
        < sum(abs(Wi(2:end))) + sum(abs(We))  );

    if triali
        triali_count = triali_count + triali;
    end

end

prob = (triale_count + triali_count)/(2*runs);

end