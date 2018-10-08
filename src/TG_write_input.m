function TG_write_input(Params, out_dir)
%myFun - Writes TG_main input to file
%
% Syntax: output = myFun(input)
%
% Long description

file = [out_dir 'results.txt'];
fileID = fopen(file, 'w');

var_names = {'N', 'f', 'tau', 'sigmae', 'sigmai', 'mue', 'mui', ...
    'frac_EV_TG', 't_min', 't_max', 't_step', 'seed'};


fprintf(fileID, datestr(now,'mm/dd/yyyy HH:MM:SS\n'));
fprintf(fileID, 'Parameter values loaded: \n');
fprintf(fileID, '---------------------------\n');
for i = 1:length(var_names)
    name = var_names{i};
    eval(['val = Params.' name ';']);
    if ischar(val)
        fprintf(fileID, '  %-12s = %s\n', name, val);
    else
        fprintf(fileID, '  %-12s = %g\n',name,val);
    end
end
fprintf(fileID, '---------------------------\n\n');


end