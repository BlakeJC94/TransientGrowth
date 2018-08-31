function TG_plot_max_growth(G_vec, t_vec, G_stats)
%TG_plot_max_growth - Plots max growth G against t
%
% Syntax: output = TG_plot_max_growth(input)
%
% Simple function for plotting the max growth as 
%   calculated from TG_get_max_growth

%% Unpack G_stats
t_opt = G_stats.t_opt;
G_max = G_stats.G_max;
G_init_slope = G_stats.G_init_slope;


%% Plot options 
fig2opt.xlabelstr = '$t$';
fig2opt.ylabelstr = ['$G(t):=\displaystyle\max_{{\bf u(t=0)}}' ...
    '\frac{{\bf u}(t)}{{\bf u}(t=0)}$'];

t_opt_str = ['$t_{opt} = $' num2str(t_opt)];
G_max_str = ['$G_{max} = $' num2str(G_max)];
G_init_slope_str = ['$\frac{dG}{dt}|_{t=0} = $' num2str(G_init_slope)];



%% Plot figure
figure;
plot(t_vec, G_vec, '--')
xlabel(fig2opt.xlabelstr,'interpreter','latex');
ylabel(fig2opt.ylabelstr,'interpreter','latex');

hold on;
% plot dashed line at peak
plot([t_opt, t_opt], [G_max, 0], 'k--');
% annotate G_stats
text(t_opt+0.1, 0.1, t_opt_str, 'Interpreter', 'latex');
text(t_opt+0.1, G_max+0.1, G_max_str, 'Interpreter', 'latex');
text(t_vec(1)+0.1, G_vec(1), G_init_slope_str, 'Interpreter', 'latex');
hold off;



end