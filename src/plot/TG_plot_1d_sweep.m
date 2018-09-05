function TG_plot_1d_sweep(sweep_results, var_str, var_vec, out_dir)
%myFun - Description
%
% Syntax: output = myFun(input)
%
% Long description

total_runs = size(sweep_results, 2);

G_max_avg = mean(sweep_results(:,:,1), 2);
G_max_std = std(sweep_results(:,:,1), 0, 2);

t_opt_avg = mean(sweep_results(:,:,2), 2);
t_opt_std = std(sweep_results(:,:,2), 0, 2);

G_init_slope_avg = mean(sweep_results(:,:,3), 2);
G_init_slope_std = std(sweep_results(:,:,3), 0, 2);



figure;
errorbar(var_vec,G_max_avg,G_max_std,'ko-'); 
xlim([min(var_vec)-0.1*range(var_vec), max(var_vec)+0.1*range(var_vec)]);
xlabel(var_str);
ylabel('G\_max')
title(['G\_max average over ' num2str(total_runs) ' runs']);
plot_export_fig(-1, [out_dir 'figures/g_max_avg'], 14, 7/5, 18);

figure;
errorbar(var_vec,t_opt_avg,t_opt_std,'ko-');
xlim([min(var_vec)-0.1*range(var_vec), max(var_vec)+0.1*range(var_vec)]);
xlabel(var_str);
ylabel('t\_opt')
title(['t\_opt average over ' num2str(total_runs) ' runs']);
plot_export_fig(-1, [out_dir 'figures/t_opt_avg'], 14, 7/5, 18);

figure;
errorbar(var_vec,G_init_slope_avg,G_init_slope_std,'ko-');        
xlim([min(var_vec)-0.1*range(var_vec), max(var_vec)+0.1*range(var_vec)]);
xlabel(var_str);
ylabel('G\_init\_slope')
title(['G\_init\_slope average over ' num2str(total_runs) ' runs']);
plot_export_fig(-1, [out_dir 'figures/G_init_slope_avg'], 14, 7/5, 18);


end