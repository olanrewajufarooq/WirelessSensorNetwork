function individual_plot_data_multi(fig_number, rounds, dims, initial_SN, sim_params, sim_name)
%PLOT_DATA_MULTI Plot of all graphs in the simulation
%   This function plots all the graphs involved in the wireless sensor
%   network. It plots the WSN and the graphs of the analytics.
%
%   INPUT PARAMETERS
%   fig_number - figure numbeer for plot.
%   rounds - the total number of rounds in the simualtion.
%   initial_SN - all sensors nodes during initiation (including routing 
%                   nodes).
%   dims - container of the dimensions of the WSN plot extremes and the
%           base station point. outputs: x_min, x_min, y_min, x_max, y_max, 
%           bs_x, bs_y
%   sim_params - container of the parameters of the data gathered after a
%                   complete imulation round. The parameters are vectors.
%                   They include: "dead nodes", "operating nodes", 
%                   "total energy", "packets", "cluster heads".
%   sim_name - the name of the wireless network simulation. 
%               Default: 'LEACH'.

if nargin < 6
    sim_name = 'LEACH';
end

colors = containers.Map( {1, 2, 3, 4, 5, 6}, {'-r', '-g', '-b', '-k', '-m', '-y'} );

i = 0;
for algorithm = 1:length(initial_SN)
    figure(fig_number + i)
    plotWSN(initial_SN(algorithm), dims, char(sim_name(algorithm)));
    i = i + 1;
end

j = 0;
label_names = ["Number of dead nodes", "Number of alive nodes", "Total energy consumed", "Number of packets received at BS", "Number of cluster heads"];

for param = ["dead nodes", "operating nodes", "total energy", "packets", "cluster heads"]
    i = i + 1;
    figure(fig_number + i)
    for algorithm=1:length(sim_params)
        sim_param = sim_params(algorithm);
        plot(1:rounds,sim_param(param),colors(algorithm),'Linewidth',2);
        hold on
    end
    xlim([0 rounds]);
    axis tight
    title( [capitalize(param), 'Per Round'] );
    xlabel 'Rounds';
    j = j + 1;
    ylabel ( capitalize(label_names(j)) );
    legend(sim_name);
    hold off
end

end

