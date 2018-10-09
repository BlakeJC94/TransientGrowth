function TG_ode_solve


% just use this to test plotting style


N = 10;
t_max = 3;



out_dir = 'output/test/wumbo/';


t_vec = linspace(0, t_max, 100);
u_mat = zeros(N, length(t_vec));

for ind = 2:length(t_vec)
    u_mat(:, ind) = t_vec(ind)^2 .* (1:N);

end


keyboard;

heatmap(1:N, t_vec, u_mat)
plot_export_fig(-1, [out_dir 'wumbojiggalo'], 14, 7/5, 18);


    
    





end