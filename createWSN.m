function [SN, rnid] = createWSN(nodes, route_nodes, dims, energy, method, seed, rn_dist)
%CREATEWSN Creation of the Wireless Sensor Network
%   This function gives the initialization of the sensor nodes, the routing
%   nodes and the base station of the wireless sensor network (WSN). It
%   also initiates the energy values and some important conditions for the
%   WSN simulation.
%
%   INPUT PARAMETERS
%   nodes - the total number of sensor and routing nodes
%   route_nodes - the number of routing nodes
%   dims - the dimensions of the WSN
%   energy - initial energy of the nodes
%   method - method of routing node formation. There are four options.
%               'equi' => equidistant from the base station, 'segcentre' =>
%               at the centre of each routing cluster, 'segedge' => at the
%               edge of each cluster, 'random' => random selection of the
%               routing node.
%   rn_dist - distance of routing node. Needed only if method = 'equi'.
%   seed - the random generation seed. Default: true. But you can pass a
%               new seed by assigning a numeric valid to the seed
%               parameter. If you don't want seeding, assign 'false'.
%
%   OUTPUT PARAMETERS
%   SN - all sensors nodes (including routing routes)
%   rnid - ids of all sensor nodes used for routing

%% Function Default Values

if route_nodes >= nodes
    error('The number of routing nodes must be less than the total number of nodes')
end

if nargin < 6
    seed = true;
    method = 'equi';
    rn_dist = dims('x_max')/4;
elseif nargin < 7
    rn_dist = 0;
    seed = true;
end

% Simulation Seed Initiation
if seed == true
    i_seed = 0;
elseif isnumeric(seed)
    i_seed = seed;
end


%% Building the sensor nodes of the WSN

SN = struct();

for i=1:nodes
        
    if seed ~= false
        rng(i_seed);
        i_seed = i_seed + 1;
    end

    SN.n(i).id = i;	% sensor's ID number
    SN.n(i).x = rand(1,1)*dims('x_max');	% X-axis coordinates of sensor node
    SN.n(i).y = rand(1,1)*dims('y_max');	% Y-axis coordinates of sensor node
    SN.n(i).E = energy;     % nodes energy levels (initially set to be equal to "ener('init')"
    SN.n(i).role = 'N';   % node acts as normal if the value is 'N', if elected as a cluster head it  gets the value 'C' (initially all nodes are normal).Nodes can also be designed as routing => 'R'
    SN.n(i).cluster = 0;	% the cluster which a node belongs to
    SN.n(i).cond = 'A';	% States the current condition of the node. when the node is operational (i.e. alive) its value = 'A' and when dead, value = 'D'
    SN.n(i).rop = 0;	% number of rounds node was operational
    SN.n(i).rleft = 0;  % rounds left for node to become available for Cluster Head election
    SN.n(i).dnc = 0;	% nodes distance from the cluster head of the cluster in which he belongs
    SN.n(i).dnb = sqrt( (dims('bs_x')-SN.n(i).x)^2 + (dims('bs_y')-SN.n(i).y)^2 );    % nodes distance from the base station
    SN.n(i).chelect = 0;	% states how many times the node was elected as a Cluster Head
    SN.n(i).rnd_chelect = 0;     % round node got elected as cluster head
    SN.n(i).chid = 0;   % node ID of the cluster head which the "i" normal node belongs to
    
end

%% Building the routing nodes of the WSN

rnid = zeros(1, route_nodes);

if strcmpi(method, 'equi')
    
    for i=1:route_nodes
        theta = i * 2*pi / route_nodes;
        route_x = rn_dist*cos(theta) + dims('bs_x');
        route_y = rn_dist*sin(theta) + dims('bs_y');
        
        % Computation of distance of each sensor node to the routing node
        % mark point
        dsr = zeros(1, nodes);
        for j = 1:nodes
            if strcmp(SN.n(j).role, 'N')             
                dist = sqrt( (route_x - SN.n(j).x)^2 + (route_y - SN.n(j).y)^2 );
                dsr(j) = dist;
            else
                dsr(j) = dims('x_max');
            end
        end
        
        % Determination of the ID of the sensor node closest to the
        % location of the desired routing node
        [~, route_id] = min( dsr(:) );
        
        SN.n(route_id).role = 'R';
        rnid(i) = route_id;
        
    end
    
elseif strcmpi(method, 'segcentre')
    
elseif strcmpi(method, 'segedge')
    
elseif strcmpi(method, 'random')
    
end

end

