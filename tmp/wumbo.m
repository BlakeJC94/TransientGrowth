function wumbo
%myFun - Description
%
% Syntax: output = myFun(input)
%
% Long description

run('parameters.m')

var_names = {'N', 'f', 'tau', 'sigmae', 'sigmai', 'mue', 'mui', ...
    'frac_EV_TG', 't_min', 't_max', 't_step', 'seed'};

% fileID = fopen('test.txt', 'w');


diary('test.txt')

fprintf('Parameter Values loaded : \n');
fprintf('---------------------------\n');
for i = 1:length(var_names)
    name = var_names{i};
    eval(['val = ' name ';']);
    fprintf('  %-12s = %g\n',name,val);
end

diary off




    
end