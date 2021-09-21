function [SN, round_params] = energy_dissipation(SN, CL_heads, round, dims, energy, rn_ids, k, round_params, method)
%ENERGY_DISSIPATION Energy dissipation function for the WSN
%   This function evaluates the energy dissipated in the sensor nodes
%   during the transmission netween the nodes to the base station of the
%   network
%
%   INPUT PARAMETERS
%   SN - all sensors nodes (including routing routes)
%   CLheads - number of cluster heads elected.
%   round - the current round in the simulation.
%   dims - container of the dimensions of the WSN plot extremes and the
%           base station point. outputs: x_min, x_min, y_min, x_max, y_max, 
%           bs_x, bs_y.
%   ener - container of the energy values needed in simulation for the
%           transceiver, amplification, aggregation. Outputs: init, tran,
%           rec, amp, agg.
%   rn_ids - ids of all sensor nodes used for routing
%   k - the number of bits transfered per packet
%   round_params - container of the parameters used to measure the
%                   performance of the simulation in a round. The params
%                   are: 'dead nodes', 'operating nodes', 'total energy', 
%                   'packets', 'stability period', 'lifetime', 
%                   'stability period round', 'lifetime round'.
%   method - the approach used in the transfer of data from normal nodes to
%               the base station. The available parameters are: 'force CH'
%               and 'shortest'. Default: 'force CH'. 'force CH' compels the
%               nodes to pass through a channel head. 'shortest' searches
%               for the minimum energy dissipation route.
%
%   OUTPUT PARAMETERS
%   SN - all sensors nodes (including routing routes)
%   round_params - container of the parameters used to measure the
%                   performance of the simulation in a round. The params
%                   are: 'dead nodes', 'operating nodes', 'total energy', 
%                   'packets', 'stability period', 'lifetime', 
%                   'stability period round', 'lifetime round'.


if nargin < 9
    method = 'force CH';
end

% Normal Nodes
for i=1:length(SN.n)
    if strcmp(SN.n(i).cond,'A') && strcmp(SN.n(i).role, 'N') && CL_heads > 0
        if SN.n(i).E > 0
            if strcmp(method, 'force CH')
            
                ETx = energy('tran')*k + energy('amp') * k * SN.n(i).dnc^2;
                SN.n(i).E = SN.n(i).E - ETx;
                round_params('total energy') = round_params('total energy') + ETx;

                % Dissipation for cluster head during reception
                if SN.n(SN.n(i).chid).E > 0 && strcmp(SN.n(SN.n(i).chid).cond, 'A') && strcmp(SN.n(SN.n(i).chid).role, 'C')
                    ERx=(energy('rec') + energy('agg'))*k;
                    round_params('total energy') = round_params('total energy') + ERx;
                    SN.n(SN.n(i).chid).E = SN.n(SN.n(i).chid).E - ERx;
                    if SN.n(SN.n(i).chid).E<=0  % if cluster heads energy depletes with reception
                        SN.n(SN.n(i).chid).cond = 'D';
                        SN.n(SN.n(i).chid).rop=round;
                        round_params('dead nodes') = round_params('dead nodes') + 1;
                        round_params('operating nodes') = round_params('operating nodes') - 1;
                    end
                
                end

            elseif strcmp(method, 'shortest')
            
            
                % Initializing the distance matrix
                dist = zeros(1, length(rn_ids) + 1);

                % Search for closest routing node or base station
                for j=1:length(rn_ids)

                    if strcmp(SN.n(rn_ids(j)).cond,'A')
                        % distance of cluster head to routing node
                        dist(j)=sqrt((SN.n(rn_ids(j)).x-SN.n(i).x)^2 + (SN.n(rn_ids(j)).y-SN.n(i).y)^2);
                    else
                        % distance of cluster head to base station
                        dist(j)=sqrt( (dims('bs_x')-SN.n(i).x)^2 + (dims('bs_y')-SN.n(i).y)^2 ) + 1;
                    end

                end

                dist(length(rn_ids) + 1) = SN.n(i).dnc;
                dist(length(rn_ids) + 2) = SN.n(i).dnb;

                [~,I]=min(dist(:)); % finds the minimum distance of node to RN

                % Direct transmission to base station
                if I == length(rn_ids) + 2 

                    ETx = (energy('tran')+energy('agg'))*k + energy('amp') * k * SN.n(i).dnb^2;
                    SN.n(i).E = SN.n(i).E - ETx;
                    round_params('total energy') = round_params('total energy') + ETx;
                    round_params('packets') = round_params('packets') + 1;

                % Transmission via cluster head
                elseif I == length(rn_ids) + 1 && SN.n(i).cluster ~= 0
                    ETx = energy('tran')*k + energy('amp') * k * SN.n(i).dnc^2;
                    SN.n(i).E = SN.n(i).E - ETx;
                    round_params('total energy') = round_params('total energy') + ETx;

                    % Dissipation for cluster head during reception
                    if SN.n(SN.n(i).chid).E > 0 && strcmp(SN.n(SN.n(i).chid).cond, 'A') && strcmp(SN.n(SN.n(i).chid).role, 'C')
                        ERx=(energy('rec') + energy('agg'))*k;
                        round_params('total energy') = round_params('total energy') + ERx;
                        SN.n(SN.n(i).chid).E = SN.n(SN.n(i).chid).E - ERx;
                        if SN.n(SN.n(i).chid).E <= 0  % if cluster heads energy depletes with reception
                            SN.n(SN.n(i).chid).cond = 'D';
                            SN.n(SN.n(i).chid).rop=round;
                            round_params('dead nodes') = round_params('dead nodes') + 1;
                            round_params('operating nodes') = round_params('operating nodes') - 1;
                        end
                    end

                % Transmission via routing
                else           

                    dnr = dist(I); % assigns the distance of node to RN
                    SN.n(i).route_id = rn_ids(I);

                    % Transmission energy to the Routing Node
                    ETx = (energy('tran')+energy('agg'))*k + energy('amp') * k * dnr^2;
                    SN.n(i).E = SN.n(i).E - ETx;
                    round_params('total energy') = round_params('total energy') + ETx;
                end  
            end
            
        end
        
        % Check for node depletion
        if SN.n(i).E<=0 % if nodes energy depletes with transmission
            round_params('dead nodes') = round_params('dead nodes') + 1;
            round_params('operating nodes') = round_params('operating nodes') - 1;
            SN.n(i).cond = 'D';
            SN.n(i).chid=0;
            SN.n(i).rop=round;
        end
        
    end        
end 

% Energy Dissipation in Cluster Head %
for i=1:length(SN.n)
    if strcmp(SN.n(i).role, 'C') &&  strcmp(SN.n(i).cond,'A')
        
        % Initializing the distance matrix
        dist = zeros(1, length(rn_ids) + 1);
        
        % Search for closest routing node or base station
        for j=1:length(rn_ids)

            if strcmp(SN.n(rn_ids(j)).cond,'A')
                % distance of cluster head to routing node
                dist(j)=sqrt((SN.n(rn_ids(j)).x-SN.n(i).x)^2 + (SN.n(rn_ids(j)).y-SN.n(i).y)^2);
            else
                % distance of cluster head to base station
                dist(j)=sqrt( (dims('bs_x')-SN.n(i).x)^2 + (dims('bs_y')-SN.n(i).y)^2 ) + 1;
            end

        end
        
        dist(length(rn_ids) + 1) = SN.n(i).dnb;

        [~,I]=min(dist(:)); % finds the minimum distance of node to RN
        
        % Direct transmission without routing
        if I == length(rn_ids) + 1 
            
            ETx = (energy('tran')+energy('agg'))*k + energy('amp') * k * SN.n(i).dnb^2;
            SN.n(i).E = SN.n(i).E - ETx;
            round_params('total energy') = round_params('total energy') + ETx;
            round_params('packets') = round_params('packets') + 1;
         
        % Transmission via routing
        else           
            
            dcr = dist(I); % assigns the distance of node to RN
            SN.n(i).route_id = rn_ids(I);
            
            % Transmission energy to the Routing Node
            ETx = (energy('tran')+energy('agg'))*k + energy('amp') * k * dcr^2;
            SN.n(i).E = SN.n(i).E - ETx;
            round_params('total energy') = round_params('total energy') + ETx;
            
            % Receiving energy at the Routing Node
            ERx=( energy('rec') + energy('agg') )*k;
            round_params('total energy') = round_params('total energy') + ERx;
            SN.n(SN.n(i).route_id).E=SN.n(SN.n(i).route_id).E - ERx;
            
            % Transmission from the routing node to the BS %
            if SN.n(SN.n(i).route_id).E>0 % Transmission to BS
                ETx= (energy('tran')+energy('agg'))*k + energy('amp') * k * SN.n(SN.n(i).route_id).dnb^2;
                SN.n(SN.n(i).route_id).E = SN.n(SN.n(i).route_id).E - ETx;
                round_params('total energy') = round_params('total energy') + ETx;
                round_params('packets') = round_params('packets') + 1;
            end
            
            if SN.n(SN.n(i).route_id).E <= 0  % if cluster heads energy depletes with reception
                SN.n(SN.n(i).route_id).cond = 'D';
                SN.n(SN.n(i).route_id).rop=round;
                round_params('dead nodes') = round_params('dead nodes') + 1;
                round_params('operating nodes') = round_params('operating nodes') - 1;
            end
            
        end
        
        if  SN.n(i).E<=0 % if cluster heads energy depletes with transmission
            round_params('dead nodes') = round_params('dead nodes') + 1;
            round_params('operating nodes') = round_params('operating nodes') - 1;
            SN.n(i).cond='D';
            SN.n(i).rop=round;
        end
        
    end
end


end

