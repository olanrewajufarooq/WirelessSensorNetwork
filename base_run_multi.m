close all;
clear;
clc;

%% Simulation Details

% Initialization Inputs
max_dimension = 100; % Maximum Dimension of the WSN Plot

initial_energy = 500e-3; % Initial node energy
transceiver_energy = 50e-9; % Energy Required for transmission and receiving of packet
ener_amp = 100e-12; % Amplification Energy
ener_agg = 100e-12; % Aggregation Energy

% Parameter Initialization
[dims, ener] = param_init(max_dimension, initial_energy, transceiver_energy, ener_agg, ener_amp);

% Simulation Parameters

n = 100; % Number of nodes
rn_old = 4; % Number of routing nodes
rn_new = 4;

rn_dist = 30;
dims('rn_dist') = rn_dist;

rounds = 1000; % Number of rounds per simulation
% sim = 3; % Number of simulations
k = 8000; % Bits transmitted per packet

% Clustering Paramters
p=0.05; % Percentage of cluster heads

%% Initialization of the WSN
[initial_SN_old, rn_ids_old] = createWSN(n, rn_old, dims, ener('init'), 'equi', 7, dims('rn_dist'));
[initial_SN_new, rn_ids_new] = createWSN(n, rn_new, dims, ener('init'), 'equi', 200, dims('rn_dist'));

%% Smiluation of the WSN
[SN_old, round_params_old, sim_params_old] = simulation_rounds(rounds, initial_SN_old, p, k, dims, ener, rn_ids_old);
[SN_new, round_params_new, sim_params_new] = simulation_rounds(rounds, initial_SN_new, p, k, dims, ener, rn_ids_new, true);

%% Graphs
plot_data(1, rounds, dims, initial_SN_old, sim_params_old)
plot_data(2, rounds, dims, initial_SN_new, sim_params_new)

%% Combined Graph
initial_SN = containers.Map( {1, 2}, {initial_SN_old, initial_SN_new} );
sim_params = containers.Map( {1, 2}, {sim_params_old, sim_params_new} );
plot_data_multi(3, rounds, dims, initial_SN, sim_params, {'LEACH', 'New ALgorithm'})

%% Individual Plots

%{
individual_plot_data(4, rounds, dims, initial_SN_old, sim_params_old, 'LEACH')
individual_plot_data(11, rounds, dims, initial_SN_new, sim_params_new, 'New Algorithm')
%}

%% Individual Plots (Multiple)

individual_plot_data_multi(4, rounds, dims, initial_SN, sim_params, {'LEACH', 'New ALgorithm'})

%% Lifetime and Stability Periods.

fprintf('\n\nLEACH Algorithm\n')
fprintf('Stability Period: %d\n', round(round_params_old('stability period'), 2))
fprintf('Stability Period Round: %d\n', round_params_old('stability period round'))
fprintf('Lifetime: %d\n', round(round_params_old('lifetime'), 2))
fprintf('Lifetime Round: %d\n', round_params_old('lifetime round'))

fprintf('\n\nNew Algorithm\n')
fprintf('Stability Period: %d\n', round(round_params_new('stability period'), 2))
fprintf('Stability Period Round: %d\n', round_params_new('stability period round'))
fprintf('Lifetime: %d\n', round(round_params_new('lifetime'), 2))
fprintf('Lifetime Round: %d\n', round_params_new('lifetime round'))