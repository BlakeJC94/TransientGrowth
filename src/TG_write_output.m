function TG_write_output(G_stats, out_dir, name, value, index)
%myFun - Writes TG_main input to file
%
% Syntax: output = myFun(input)
%
% Long description

file = [out_dir 'results.txt'];
fileID = fopen(file, 'a');

if nargin > 2
    fprintf(fileID, 'Results (%s = %g, %g): \n', name, value, index);
else
    fprintf(fileID, 'Results : \n');
end

fprintf(fileID, '---------------------------\n');
fprintf(fileID,'  %-12s = %g\n','G_max', G_stats.G_max);
fprintf(fileID,'  %-12s = %g\n','t_opt', G_stats.t_opt);
fprintf(fileID,'  %-12s = %g\n','G_init_slope', G_stats.G_init_slope);
fprintf(fileID, '---------------------------\n\n');




end