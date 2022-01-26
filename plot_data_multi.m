function plot_data_multi(fig_number, rounds, dims, initial_SN, sim_params, sim_name)
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

plot_num = ceil( 5 + length(initial_SN) );
row_num = ceil(plot_num / 3);

i = 1;
figure(fig_number)
for algorithm = 1:length(initial_SN)
    subplot(row_num, 3, i)
    plotWSN(initial_SN(algorithm), dims, char(sim_name(algorithm)));
    i = i + 1;
end


for param = ["dead nodes", "operating nodes", "total energy", "packets", "cluster heads"]
    i = i + 1;
    subplot(row_num, 3, i)
    
    for algorithm=1:length(sim_params)
        sim_param = sim_params(algorithm);
        plot(1:rounds,sim_param(param),colors(algorithm),'Linewidth',2);
        hold on
    end
    
    xlim([0 rounds]);
    axis tight
    title( [capitalize(param), 'Per Round'] );
    xlabel 'Rounds';
    
    ylabel ( capitalize(param) );
    legend(sim_name);
end
hold off

end

