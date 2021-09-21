function plot_data(fig_number, rounds, dims, initial_SN, sim_params, sim_name)
%PLOT_DATA Plot of all graphs in the simulation
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

i = 1;
figure(fig_number)
subplot(2, 3, i)
plotWSN(initial_SN, dims, sim_name)


for param = ["dead nodes", "operating nodes", "total energy", "packets", "cluster heads"]
    i = i + 1;
    subplot(2, 3, i)
    plot(1:rounds,sim_params(param),'-r','Linewidth',2);
    xlim([0 rounds]);
    axis tight
    title( [capitalize(param), 'Per Round'] );
    xlabel 'Rounds';
    ylabel ( capitalize(param) );
    legend(sim_name);
end


end

