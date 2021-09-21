function [SN, CL, CLheads] = cluster_head_election(SN, p, round, dims, seed)
%CLUSTER_HEAD_ELECTION Election of the Cluster Heads
%   This function gives the selection of the cluster heads in a wireless
%   sensor network (WSN) per round. 
%
%   INPUT PARAMETERS
%   SN - all sensors nodes (including routing routes)
%   p - the percentage of cluster heads desired
%   round - the round number of the simulation
%   dims - the dimensions of the WSN
%   seed - the random generation seed. Default: true. But you can pass a
%               new seed by assigning a numeric valid to the seed
%               parameter. If you don't want seeding, assign 'false'.
%
%   OUTPUT PARAMETERS
%   SN - all sensors nodes (including routing routes)
%   CL - all cluster nodes
%   CLheads - number of cluster heads elected.

if nargin < 5
    seed = true;
end

% Election Seed
if seed == true
    i_seed = 0;
elseif isnumeric(seed)
    i_seed = seed;
end


% Threshold for cluster election
t = p/( 1 - p*(mod(round, 1/p)) );
% Re-election Value
tleft = mod(round, 1/p);

% Initiating the cluster heads
CLheads = 0; % Initializing the cluster head count
CL = struct(); % The struct containing the cluster head data

for i=1:length(SN.n)
    
    if seed ~= false
        rng(i_seed);
        i_seed = i_seed + 1;
    end
    
    if SN.n(i).role ~= 'R'
    
        SN.n(i).cluster=0;    % reseting cluster in which the node belongs to
        SN.n(i).role='N';       % reseting node role
        SN.n(i).chid=0;       % reseting cluster head id

        % If node has been initially elected as cluster head
        if SN.n(i).rleft > 0
           SN.n(i).rleft = SN.n(i).rleft-1;
        end

        % If Node is eligible to be elected as channel head
        if (SN.n(i).E > 0) && (SN.n(i).rleft == 0)
            generate = rand;

            if generate < t
                SN.n(i).role = 'C';	% assigns the node role of a cluster head
                SN.n(i).rnd_chelect = round;	% Assigns the round that the cluster head was elected to the data table
                SN.n(i).chelect = SN.n(i).chelect + 1;   
                SN.n(i).rleft = 1/p - tleft;    % rounds for which the node will be unable to become a CH
                SN.n(i).dcb = sqrt( (dims('bs_x')-SN.n(i).x)^2 + (dims('bs_y')-SN.n(i).y)^2 ); % calculates the distance between the base station and the cluster head
                CLheads = CLheads + 1;
                SN.n(i).cluster = CLheads; % cluster of which the node got elected to be cluster head
                CL.n(CLheads).x=SN.n(i).x; % X-axis coordinates of elected cluster head
                CL.n(CLheads).y=SN.n(i).y; % Y-axis coordinates of elected cluster head
                CL.n(CLheads).id=i; % Assigns the node ID of the newly elected cluster head to an array
            end
            
        end
        
    end
end

end

