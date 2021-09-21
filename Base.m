close all;
clear;
clc;

%% Network Establishment Parameters 

% Number of simulations
sim_n=3;

%%% Area of Operation 

% Field Dimensions in meters %
xm=100;
ym=100;

x=0; % added for better display results of the plot
y=0; % added for better display results of the plot

% Coordinates of the Sink (location is predetermined in this simulation) %
sinkx=xm/2;
sinky=ym/2;

%%% Energy Values

% Initial Energy of a Node (in Joules) % 
Eo=500*10^(-3); % units in Joules

% Energy required to run circuity (both for transmitter and receiver) %
Eelec=50*10^(-9); % units in Joules/bit
ETx=50*10^(-9); % units in Joules/bit
ERx=50*10^(-9); % units in Joules/bit

% Transmit Amplifier Types %
Eamp=100*10^(-12); % units in Joules/bit/m^2 (amount of energy spent by the amplifier to transmit the bits)

% Data Aggregation Energy %
EDA=5*10^(-9); % units in Joules/bit

% Size of data package %
k=800; % units in bits

% Suggested percentage of cluster head %
p=0.05; % a 5 percent of the total amount of nodes used in the network is proposed to give good results

% Maximum Number of Rounds %
max_rnd = 6000;


%% Parameters for LEACH Algorithm

% Number of Nodes in the field %
n=104;

% Plot Details for the Routing Circle
drnbs = 30; % distance of routing node to base station
theta = linspace(0,2*pi);
xcir = drnbs*cos(theta) + sinkx;
ycir = drnbs*sin(theta) + sinky;
routemark = [0 pi/2 pi 3*pi/2];
xcirmark = drnbs*cos(routemark) + sinkx;
ycirmark = drnbs*sin(routemark) + sinky;


%% Parameters for the New Algorithm 

% Number of Nodes in the field %
new_n=104;
% Number of routers %
nr = 8;

%% Creation of the Wireless Sensor Network

% Initialize random seed
seed = 0;

for sim_k=1:sim_n
    rng(seed)
    
    close all;
    
    %% General Initiazations
    % Round of Operation %
    rnd = 0;
    % Moving Average method
    s = 0.1; % smoothing value
    
    
    %% Initialization for LEACH Algorithm
    
    % Current Number of operating Nodes %
    operating_nodes = n;
    % Number of Dead Nodes in the beggining %
    dead_nodes=0;
    % Number of packets received
    pkt = 0;
    % Total energy consumed
    tot_energy = 0;
    % Average energy consumed
    avg_energy = 0;
    % Moving Average method
    mvg_avg_energy = 0;
    
    
    %% Initialization for the New ALgorithm
    
    % Number of Dead Nodes in the beggining %
    new_dead_nodes = 0;
    % Current Number of operating Nodes %
    new_operating_nodes=new_n;
    % Number of packets received
    new_pkt = 0;
    % Total energy consumed
    new_tot_energy = 0;
    % Average energy consumed
    new_avg_energy = 0;
    % Moving Average method
    new_mvg_avg_energy = 0;
    
    
    %% Plotting the WSN for LEACH Algorithm
    for i=1:n
        
        %Calling random seed
        seed = seed + 1
        rng(seed)
        
        SN(sim_k,i).id=i;	% sensor's ID number
        SN(sim_k,i).x=rand(1,1)*xm;	% X-axis coordinates of sensor node
        SN(sim_k,i).y=rand(1,1)*ym;	% Y-axis coordinates of sensor node
        SN(sim_k,i).E=Eo;     % nodes energy levels (initially set to be equal to "Eo"
        SN(sim_k,i).role=0;   % node acts as normal if the value is '0', if elected as a cluster head it  gets the value '1' (initially all nodes are normal)
        SN(sim_k,i).cluster=0;	% the cluster which a node belongs to
        SN(sim_k,i).cond=1;	% States the current condition of the node. when the node is operational its value is =1 and when dead =0
        SN(sim_k,i).rop=0;	% number of rounds node was operational
        SN(sim_k,i).rleft=0;  % rounds left for node to become available for Cluster Head election
        SN(sim_k,i).dtch=0;	% nodes distance from the cluster head of the cluster in which he belongs
        SN(sim_k,i).dts=0;    % nodes distance from the sink
        SN(sim_k,i).tel=0;	% states how many times the node was elected as a Cluster Head
        SN(sim_k,i).rn=0;     % round node got elected as cluster head
        SN(sim_k,i).chid=0;   % node ID of the cluster head which the "i" normal node belongs to


        hold on;
        figure(1)
        subplot(1,2,1)
        plot(x,y,xm,ym,SN(sim_k,i).x,SN(sim_k,i).y,'om');
        plot(xcir,ycir,'-b','LineWidth',3);
        plot([0 100],[50 50],'-r','LineWidth',2);
        plot([50 50],[0 100],'-r','LineWidth',2);
        plot(xcirmark,ycirmark,'dr','MarkerFaceColor','w','MarkerEdgeColor','r','MarkerSize',10);
        plot(sinkx,sinky,'*k','LineWidth',6,'MarkerSize',15);
        title ({'LEACH'; 'Wireless Sensor Network';})
        xlabel '(m)';
        ylabel '(m)';

    end

    % Setting Up Routing Node %
    new_drnbs = drnbs;
    for i=(n+1):(n+4)
        theta = pi/4 + (i-n) * 2* pi/4;

        SN(sim_k,i).id = i;
        SN(sim_k,i).x=new_drnbs*cos(theta) + sinkx;	% X-axis coordinates of sensor node
        SN(sim_k,i).y=new_drnbs*sin(theta) + sinky;	% Y-axis coordinates of sensor node
        SN(sim_k,i).E=Eo;     % nodes energy levels (initially set to be equal to "Eo"
        SN(sim_k,i).cond=1;	% States the current condition of the node. when the node is operational its value is =1 and when dead =0
        SN(sim_k,i).dts=new_drnbs;    % routing nodes distance from the sink
    end

    %%% Plotting the WSN for New Algorithm
    for i=1:new_n
        
        %Calling random seed
        seed = seed + 1
        rng(seed)
                
        nSN(sim_k,i).id=i;	% sensor's ID number
        nSN(sim_k,i).x=rand(1,1)*xm;	% X-axis coordinates of sensor node
        nSN(sim_k,i).y=rand(1,1)*ym;	% Y-axis coordinates of sensor node
        nSN(sim_k,i).E=Eo;     % nodes energy levels (initially set to be equal to "Eo"
        nSN(sim_k,i).role=0;   % node acts as normal if the value is '0', if elected as a cluster head it  gets the value '1' (initially all nodes are normal)
        nSN(sim_k,i).cluster=0;	% the cluster which a node belongs to
        nSN(sim_k,i).cond=1;	% States the current condition of the node. when the node is operational its value is =1 and when dead =0
        nSN(sim_k,i).rop=0;	% number of rounds node was operational
        nSN(sim_k,i).rleft=0;  % rounds left for node to become available for Cluster Head election
        nSN(sim_k,i).dtch=0;	% nodes distance from the cluster head of the cluster in which he belongs
        nSN(sim_k,i).dts=0;    % nodes distance from the sink
        nSN(sim_k,i).tel=0;	% states how many times the node was elected as a Cluster Head
        nSN(sim_k,i).rn=0;     % round node got elected as cluster head
        nSN(sim_k,i).chid=0;   % node ID of the cluster head which the "i" normal node belongs to


        hold on;
        subplot(1,2,2)
        plot(x,y,xm,ym,nSN(sim_k,i).x,nSN(sim_k,i).y,'om');
        plot([0 100],[sinkx+15, sinkx+15],'-r','LineWidth',3);
        plot([0 100],[sinkx-15, sinkx-15],'-r','LineWidth',3);
        plot([sinky+15, sinky+15],[0 100],'-r','LineWidth',3);
        plot([sinky-15, sinky-15],[0 100],'-r','LineWidth',3);
        plot(sinkx,sinky,'*k','LineWidth',6,'MarkerSize',15);
        title ({'New Algorithm'; 'Wireless Sensor Network';})
        xlabel '(m)';
        ylabel '(m)';

    end

    % Setting Up Routing Node %
    new_drnbs = 25;
    for i=(new_n+1):(new_n+nr)
        theta = (i-new_n) * 2*pi/(nr);

        nSN(sim_k,i).id = i;
        nSN(sim_k,i).x=new_drnbs*cos(theta) + sinkx;	% X-axis coordinates of sensor node
        nSN(sim_k,i).y=new_drnbs*sin(theta) + sinky;	% Y-axis coordinates of sensor node
        nSN(sim_k,i).E=Eo;     % nodes energy levels (initially set to be equal to "Eo"
        nSN(sim_k,i).cond=1;	% States the current condition of the node. when the node is operational its value is =1 and when dead =0
        nSN(sim_k,i).dts=new_drnbs;    % routing nodes distance from the sink

        hold on;
        plot(nSN(sim_k,i).x,nSN(sim_k,i).y,'dr','MarkerFaceColor','w','MarkerEdgeColor','b','MarkerSize',10);

    end



    %% Phase Initial Set-Up

    tic
    sp_check = 1; % Stability Period Check On (LEACH)
    lt_check = 1; % Lifetime Check On (LEACH)
    new_sp_check = 1; % Stability Period Check On (New Algorithm)
    new_lt_check = 1; % Lifetime Check On (New Algorithm)

    % Runing nth round
    while rnd <= max_rnd

        % Displays Current Round %     
        fprintf('simulation = %d : round = %d \n', sim_k, rnd);     

        % Threshold Value %
        t=p/(1-p*(mod(rnd,1/p)));
        % Re-election Value %
        tleft=mod(rnd,1/p);

        % Reseting Previous Amount Of Cluster Heads In the Network %
        CLheads=0;
        % Reseting Previous Amount Of Energy Consumed In the Network on the Previous Round %
        energy=0;
        % Reseting Previous Amount Of Cluster Heads In the Network %
        new_CLheads=0;
        % Reseting Previous Amount Of Energy Consumed In the Network on the Previous Round %
        new_energy=0;




        % Cluster Heads Election %



        % LEACH %

        for i=1:n
            SN(sim_k,i).cluster=0;    % reseting cluster in which the node belongs to
            SN(sim_k,i).role=0;       % reseting node role
            SN(sim_k,i).chid=0;       % reseting cluster head id

            % If node has been initially elected as cluster head
            if SN(sim_k,i).rleft>0
               SN(sim_k,i).rleft=SN(sim_k,i).rleft-1;
            end

            % If Node is eligible to be elected as channel head
            if (SN(sim_k,i).E>0) && (SN(sim_k,i).rleft==0)
                generate = rand;

                if generate< t
                SN(sim_k,i).role = 1;	% assigns the node role of a cluster head
                SN(sim_k,i).rn = rnd;	% Assigns the round that the cluster head was elected to the data table
                SN(sim_k,i).tel = SN(sim_k,i).tel + 1;   
                SN(sim_k,i).rleft = 1/p-tleft;    % rounds for which the node will be unable to become a CH
                SN(sim_k,i).dts=sqrt((sinkx-SN(sim_k,i).x)^2 + (sinky-SN(sim_k,i).y)^2); % calculates the distance between the sink and the cluster hea
                CLheads=CLheads+1;	% sum of cluster heads that have been elected 
                SN(sim_k,i).cluster = CLheads; % cluster of which the node got elected to be cluster head
                CL(CLheads).x=SN(sim_k,i).x; % X-axis coordinates of elected cluster head
                CL(CLheads).y=SN(sim_k,i).y; % Y-axis coordinates of elected cluster head
                CL(CLheads).id=i; % Assigns the node ID of the newly elected cluster head to an array
                end

            end
        end  
        % Fixing the size of "CL" array %
        CL=CL(1:CLheads);



        % The New Algorithm %

        for i=1:new_n
            nSN(sim_k,i).cluster=0;    % reseting cluster in which the node belongs to
            nSN(sim_k,i).role=0;       % reseting node role
            nSN(sim_k,i).chid=0;       % reseting cluster head id

            % If node has been initially elected as cluster head
            if nSN(sim_k,i).rleft>0
               nSN(sim_k,i).rleft=nSN(sim_k,i).rleft-1;
            end

            % If Node is eligible to be elected as channel head
            if (nSN(sim_k,i).E>0) && (nSN(sim_k,i).rleft==0)
                generate = rand;
                if generate< t
                nSN(sim_k,i).role = 1;	% assigns the node role of a cluster head
                nSN(sim_k,i).rn = rnd;	% Assigns the round that the cluster head was elected to the data table
                nSN(sim_k,i).tel = nSN(sim_k,i).tel + 1;   
                nSN(sim_k,i).rleft = 1/p-tleft;    % rounds for which the node will be unable to become a CH
                nSN(sim_k,i).dts=sqrt((sinkx-nSN(sim_k,i).x)^2 + (sinky-nSN(sim_k,i).y)^2); % calculates the distance between the sink and the cluster hea
                new_CLheads=new_CLheads+1;	% sum of cluster heads that have been elected 
                nSN(sim_k,i).cluster = new_CLheads; % cluster of which the node got elected to be cluster head
                nCL(new_CLheads).x=nSN(sim_k,i).x; % X-axis coordinates of elected cluster head
                nCL(new_CLheads).y=nSN(sim_k,i).y; % Y-axis coordinates of elected cluster head
                nCL(new_CLheads).id=i; % Assigns the node ID of the newly elected cluster head to an array
                end
            end
        end
        % Fixing the size of "CL" array %
        nCL=nCL(1:new_CLheads);



        % Grouping the Nodes into Clusters & calculating the distance between node and cluster head %


        % LEACH

        for i=1:n
            if  (SN(sim_k,i).role==0) && (SN(sim_k,i).E>0) && (CLheads>0) % if node is normal
                for m=1:CLheads
                d(m)=sqrt((CL(m).x-SN(sim_k,i).x)^2 + (CL(m).y-SN(sim_k,i).y)^2);
                % we calculate the distance 'd' between the sensor node that is
                % transmitting and the cluster head that is receiving with the following equation+ 
                % d=sqrt((x2-x1)^2 + (y2-y1)^2) where x2 and y2 the coordinates of
                % the cluster head and x1 and y1 the coordinates of the transmitting node
                end
            d=d(1:CLheads); % fixing the size of "d" array
            [M,I]=min(d(:)); % finds the minimum distance of node to CH
            [Row, Col] = ind2sub(size(d),I); % displays the Cluster Number in which this node belongs too
            SN(sim_k,i).cluster = Col; % assigns node to the cluster
            SN(sim_k,i).dtch = d(Col); % assigns the distance of node to CH
            SN(sim_k,i).chid = CL(Col).id;
            end
        end


        % New Algorithm %
        for i=1:new_n
            if  (nSN(sim_k,i).role==0) && (nSN(sim_k,i).E>0) && (new_CLheads>0) % if node is normal
                for m=1:new_CLheads
                new_d(m)=sqrt((nCL(m).x-nSN(sim_k,i).x)^2 + (nCL(m).y-nSN(sim_k,i).y)^2);
                % we calculate the distance 'd' between the sensor node that is
                % transmitting and the cluster head that is receiving with the following equation+ 
                % d=sqrt((x2-x1)^2 + (y2-y1)^2) where x2 and y2 the coordinates of
                % the cluster head and x1 and y1 the coordinates of the transmitting node
                end
            new_d=new_d(1:new_CLheads); % fixing the size of "d" array
            [M,I]=min(new_d(:)); % finds the minimum distance of node to CH
            [Row, Col] = ind2sub(size(d),I); % displays the Cluster Number in which this node belongs too
            nSN(sim_k,i).cluster = Col; % assigns node to the cluster
            nSN(sim_k,i).dtch = new_d(Col); % assigns the distance of node to CH
            nSN(sim_k,i).chid = nCL(Col).id;
            end
        end



        %%%%%% Energy Dissipation for normal nodes %%%%%%


        % LEACH %

        for i=1:n
           if (SN(sim_k,i).cond==1) && (SN(sim_k,i).role==0) && (CLheads>0)
            if SN(sim_k,i).E>0
                ETx= Eelec*k + Eamp * k * SN(sim_k,i).dtch^2;
                SN(sim_k,i).E=SN(sim_k,i).E - ETx;
                energy=energy+ETx;

            % Dissipation for cluster head during reception
            if SN(SN(sim_k,i).chid).E>0 && SN(SN(sim_k,i).chid).cond==1 && SN(SN(sim_k,i).chid).role==1
                ERx=(Eelec+EDA)*k;
                energy=energy+ERx;
                SN(SN(sim_k,i).chid).E=SN(SN(sim_k,i).chid).E - ERx;
                 if SN(SN(sim_k,i).chid).E<=0  % if cluster heads energy depletes with reception
                    SN(SN(sim_k,i).chid).cond=0;
                    SN(SN(sim_k,i).chid).rop=rnd;
                    dead_nodes=dead_nodes +1;
                    operating_nodes=operating_nodes -1;
                 end
            end
            end
            if SN(sim_k,i).E<=0       % if nodes energy depletes with transmission
            dead_nodes = dead_nodes +1;
            operating_nodes = operating_nodes - 1;
            SN(sim_k,i).cond=0;
            SN(sim_k,i).chid=0;
            SN(sim_k,i).rop=rnd;
            end
          end
        end 

        % The New ALgorithm #

        for i=1:new_n
           if (nSN(sim_k,i).cond==1) && (nSN(sim_k,i).role==0) && (new_CLheads>0)
            if nSN(sim_k,i).E>0
                ETx= Eelec*k + Eamp * k * nSN(sim_k,i).dtch^2;
                nSN(sim_k,i).E=nSN(sim_k,i).E - ETx;
                new_energy=new_energy+ETx;

            % Dissipation for cluster head during reception
            if nSN(nSN(sim_k,i).chid).E>0 && nSN(nSN(sim_k,i).chid).cond==1 && nSN(nSN(sim_k,i).chid).role==1
                ERx=(Eelec+EDA)*k;
                new_energy=new_energy+ERx;
                nSN(nSN(sim_k,i).chid).E=nSN(nSN(sim_k,i).chid).E - ERx;
                 if nSN(nSN(sim_k,i).chid).E<=0  % if cluster heads new_energy depletes with reception
                    nSN(nSN(sim_k,i).chid).cond=0;
                    nSN(nSN(sim_k,i).chid).rop=rnd;
                    new_dead_nodes=new_dead_nodes +1;
                    new_operating_nodes=new_operating_nodes -1;
                 end
            end
           end

            if nSN(sim_k,i).E<=0       % if nodes new_energy depletes with transmission
            new_dead_nodes = new_dead_nodes +1;
            new_operating_nodes = new_operating_nodes - 1;
            nSN(sim_k,i).cond=0;
            nSN(sim_k,i).chid=0;
            nSN(sim_k,i).rop=rnd;
            end
          end
        end



       % Energy Dissipation for cluster head nodes %



       % LEACH %

       for i=1:n
         if (SN(sim_k,i).cond==1)  && (SN(sim_k,i).role==1)
             if SN(sim_k,i).E>0 && SN(sim_k,i).dts<=r % Direct transmission without routing
                ETx= (Eelec+EDA)*k + Eamp * k * SN(sim_k,i).dts^2;
                SN(sim_k,i).E=SN(sim_k,i).E - ETx;
                energy=energy+ETx;
                pkt = pkt+1;

             elseif SN(sim_k,i).E>0 && SN(sim_k,i).dts>r % Transmission via routing nodes
                 for m=(n+1):(n+4)
                     if SN(m).cond==1
                        d(m)=sqrt((SN(m).x-SN(sim_k,i).x)^2 + (SN(m).y-SN(sim_k,i).y)^2);
                     else
                         d(m)=sqrt((sinkx-SN(sim_k,i).x)^2 + (sinky-SN(sim_k,i).y)^2);
                     end
                    % we calculate the distance 'd' between the sensor node that is
                    % transmitting and the routing node that is receiving with the following equation+ 
                    % d=sqrt((x2-x1)^2 + (y2-y1)^2) where x2 and y2 the coordinates of
                    % the routing node and x1 and y1 the coordinates of the transmitting node
                    % if the routing node is dead, it transmits to the base
                    % station
                 end
                d=d(1:4); % fixing the size of "d" array
                [M,I]=min(d(:)); % finds the minimum distance of node to RN
                [Row, Col] = ind2sub(size(d),I); % displays the Router Number in which this node belongs too
                dtrn = d(Col); % assigns the distance of node to RN
                SN(sim_k,i).rnid = n+Col;
                % Transmission energy to the Routing Node
                ETx= (Eelec+EDA)*k + Eamp * k * dtrn^2;
                SN(sim_k,i).E=SN(sim_k,i).E - ETx;
                energy=energy+ETx;
                % Receiving energy at the Routing Node
                ERx=(Eelec+EDA)*k;
                energy=energy+ERx;
                SN(SN(sim_k,i).rnid).E=SN(SN(sim_k,i).rnid).E - ERx;
                % Transmission from the routing node to the BS %
                if SN(SN(sim_k,i).rnid).E>0 % Transmission to BS
                    ETx= (Eelec+EDA)*k + Eamp * k * SN(SN(sim_k,i).rnid).dts^2;
                    SN(SN(sim_k,i).rnid).E=SN(SN(sim_k,i).rnid).E - ETx;
                    energy=energy+ETx;
                    pkt = pkt+1;
                end
                 if SN(SN(sim_k,i).rnid).E<=0  % if cluster heads energy depletes with reception
                    SN(SN(sim_k,i).rnid).cond=0;
                 end
             end
             if  SN(sim_k,i).E<=0 % if cluster heads energy depletes with transmission
                 dead_nodes=dead_nodes +1;
                 operating_nodes = operating_nodes -1;
                 SN(sim_k,i).cond=0;
                 SN(sim_k,i).rop=rnd;
             end
         end
       end



       % The New Algorithm %
       for i=1:new_n
         if (nSN(sim_k,i).cond==1)  && (nSN(sim_k,i).role==1)
             if nSN(sim_k,i).E>0 && nSN(sim_k,i).dts<=r % Direct transmission without routing
                ETx= (Eelec+EDA)*k + Eamp * k * nSN(sim_k,i).dts^2;
                nSN(sim_k,i).E=nSN(sim_k,i).E - ETx;
                new_energy=new_energy+ETx;
                new_pkt = new_pkt+1;
             elseif nSN(sim_k,i).E>0 && nSN(sim_k,i).dts>r % Transmission via routing nodes
                 for m=(new_n+1):(new_n+nr)
                     if nSN(m).cond==1
                        d(m)=sqrt((nSN(m).x-nSN(sim_k,i).x)^2 + (nSN(m).y-nSN(sim_k,i).y)^2);
                     else
                         d(m)=sqrt((sinkx-nSN(sim_k,i).x)^2 + (sinky-nSN(sim_k,i).y)^2);
                     end
                    % we calculate the distance 'd' between the sensor node that is
                    % transmitting and the routing node that is receiving with the following equation+ 
                    % d=sqrt((x2-x1)^2 + (y2-y1)^2) where x2 and y2 the coordinates of
                    % the routing node and x1 and y1 the coordinates of the transmitting node
                    % if the routing node is dead, it transmits to the base
                    % station
                 end
                d=d(1:nr); % fixing the size of "d" array
                [M,I]=min(d(:)); % finds the minimum distance of node to RN
                [Row, Col] = ind2sub(size(d),I); % displays the Router Number in which this node belongs too
                dtrn = d(Col); % assigns the distance of node to RN
                nSN(sim_k,i).rnid = new_n+Col;
                % Transmission new_energy to the Routing Node
                ETx= (Eelec+EDA)*k + Eamp * k * dtrn^2;
                nSN(sim_k,i).E=nSN(sim_k,i).E - ETx;
                new_energy=new_energy+ETx;
                % Receiving new_energy at the Routing Node
                ERx=(Eelec+EDA)*k;
                new_energy=new_energy+ERx;
                nSN(nSN(sim_k,i).rnid).E=nSN(nSN(sim_k,i).rnid).E - ERx; 
                % Transmission from the routing node to the BS
                if nSN(nSN(sim_k,i).rnid).E>0 % Transmission to BS
                    ETx= (Eelec+EDA)*k + Eamp * k * nSN(nSN(sim_k,i).rnid).dts^2;
                    nSN(nSN(sim_k,i).rnid).E=nSN(nSN(sim_k,i).rnid).E - ETx;
                    new_energy=new_energy+ETx;
                    new_pkt = new_pkt+1;
                end
                 if nSN(nSN(sim_k,i).rnid).E<=0  % if cluster heads new_energy depletes with reception
                    nSN(nSN(sim_k,i).rnid).cond=0;
                 end
             end
             if  nSN(sim_k,i).E<=0 % if cluster heads new_energy depletes with transmission
                 new_dead_nodes=new_dead_nodes +1;
                 new_operating_nodes = new_operating_nodes -1;
                 nSN(sim_k,i).cond=0;
                 nSN(sim_k,i).rop=rnd;
             end
         end
       end
       

       % Next Round %
       rnd=rnd +1;
       
       % Total and Average Energy
       tot_energy = tot_energy + energy;
       new_tot_energy = new_tot_energy + new_energy;

       avg_energy = tot_energy/n;
       new_avg_energy = new_tot_energy/new_n;
       
       % Residual Energy
       res_energy = 0;
       for i=1:n
           res_energy = res_energy + SN(sim_k,i).E;
       end
       avg_res_energy = res_energy / n;

       new_res_energy = 0;
       for i=1:new_n
           new_res_energy = new_res_energy + nSN(sim_k,i).E;
       end
       new_avg_res_energy = new_res_energy / new_n;
       
       %%%%%% Collating Data %%%%%%%%


       % LEACH %

       op(sim_k, rnd)=operating_nodes;
       dn(sim_k, rnd)=dead_nodes;
       CHnum(sim_k, rnd)=CLheads;
       packet(sim_k, rnd)=pkt;
       total_energy(sim_k, rnd)=tot_energy;
       average_energy(sim_k, rnd)=avg_energy;
       mvg_average_energy(sim_k, rnd)=mvg_avg_energy;
       average_res_energy(sim_k, rnd)=avg_res_energy;


       % New ALgorithm %

       new_op(sim_k, rnd)=new_operating_nodes;
       new_dn(sim_k, rnd)=new_dead_nodes;
       new_CHnum(sim_k, rnd)=new_CLheads;
       new_packet(sim_k, rnd)=new_pkt;
       new_total_energy(sim_k, rnd)=new_tot_energy;
       new_average_energy(sim_k, rnd)=new_avg_energy;
       new_mvg_average_energy(sim_k, rnd)=new_mvg_avg_energy;
       new_average_res_energy(sim_k, rnd)=new_avg_res_energy;

       % CHECKS %

       if operating_nodes < n && sp_check == 1
           stability_period = toc;
           sp_rnd=rnd;
           sp_check = 0;
       end

       if operating_nodes == 0 && lt_check == 1
           lifetime = toc;
           lt_rnd=rnd;
           lt_check = 0;
       end

       if rnd == max_rnd && sp_check == 1
           stability_period = toc;
           sp_rnd=rnd;
       end

       if rnd == max_rnd && lt_check == 1
           lifetime = toc;
           lt_rnd=rnd;
       end

       if new_operating_nodes < n && new_sp_check == 1
           new_stability_period = toc;
           new_sp_rnd=rnd;
           new_sp_check = 0;
       end

       if new_operating_nodes == 0 && new_lt_check == 1
           new_lifetime = toc;
           new_lt_rnd=rnd;
           new_lt_check = 0;
       end

       if rnd == max_rnd && new_sp_check == 1
           new_stability_period = toc;
           new_sp_rnd=rnd;
       end

       if rnd == max_rnd && new_lt_check == 1
           new_lifetime = toc;
           new_lt_rnd=rnd;
       end
    end
    
    sim_lifetime(sim_k) = lifetime;
    sim_stability_period(sim_k) = stability_period;
    sim_sp_rnd(sim_k) = sp_rnd;
    sim_lt_rnd(sim_k) = lt_rnd;
    
    sim_new_lifetime(sim_k) = new_lifetime;
    sim_new_stability_period(sim_k) = new_stability_period;
    sim_new_sp_rnd(sim_k) = new_sp_rnd;
    sim_new_lt_rnd(sim_k) = new_lt_rnd;
    
end


%% Data Analysis

mean_lifetime = mean(sim_lifetime);
mean_new_lifetime = mean(sim_new_lifetime);
fprintf('LIFETIME \n LEACH = %.2f \n New Algorithm = %.2f \n', mean_lifetime, mean_new_lifetime);

mean_stability_period = mean(sim_stability_period);
mean_new_stability_period = mean(sim_new_stability_period);
fprintf('STABILITY PERIOD \n LEACH = %.2f \n New Algorithm = %.2f \n', mean_stability_period, mean_new_stability_period);

mean_sp_rnd = round(mean(sim_sp_rnd), 0);
mean_new_sp_rnd = round(mean(sim_new_sp_rnd), 0);
fprintf('STABILITY PERIOD ROUND \n LEACH = %d \n New Algorithm = %d \n', mean_sp_rnd, mean_new_sp_rnd);

mean_lt_rnd = round(mean(sim_lt_rnd), 0);
mean_new_lt_rnd = round(mean(sim_new_lt_rnd), 0);
fprintf('LIFETIME ROUND \n LEACH = %d \n New Algorithm = %d \n', mean_lt_rnd, mean_new_lt_rnd);

% LEACH %

sim_op = mean(op);
sim_dn = mean(dn);
sim_CHnum = mean(CHnum);
sim_packet = mean(packet);
sim_total_energy = mean(total_energy);
sim_average_energy = mean(average_energy);
sim_average_res_energy = mean(average_res_energy);

% New ALgorithm %

sim_new_op = mean(new_op);
sim_new_dn = mean(new_dn);
sim_new_CHnum = mean(new_CHnum);
sim_new_packet = mean(new_packet);
sim_new_total_energy = mean(new_total_energy);
sim_new_average_energy = mean(new_average_energy);
sim_new_average_res_energy = mean(new_average_res_energy);

%% GRAPHS

if n > new_n
    ymax = n;
else
    ymax = new_n;
end


%%% Plotting Simulation Result "Dead Nodes Per Round"
figure(2)
plot(1:rnd,sim_dn(1:rnd),'-r','Linewidth',2);
hold on
plot(1:rnd,sim_new_dn(1:rnd),'-b','Linewidth',2);
xlim([0 max_rnd]);
ylim([0 ymax]);
title 'Dead Nodes per Round';
xlabel 'Rounds';
ylabel 'Dead Nodes';
hold off
legend('LEACH','New Algorithm');

%%% Plotting Simulation Results "Operating Nodes per Round"
figure(3)
plot(1:rnd,sim_op(1:rnd),'-r','Linewidth',2);
hold on
plot(1:rnd,sim_new_op(1:rnd),'-b','Linewidth',2);
xlim([0 max_rnd]);
ylim([0 ymax]);
title 'Operating Nodes per Round';
xlabel 'Rounds';
ylabel 'Live Nodes';
hold off
legend('LEACH','New Algorithm');

%%% Plotting Simulation Results "Cluster Heads per Round" 
figure(4)
plot(1:rnd,sim_CHnum(1:rnd),'-r','Linewidth',0.1);
hold on
plot(1:rnd,sim_new_CHnum(1:rnd),'-b','Linewidth',0.1);
xlim([0 max_rnd]);
title 'Cluster Heads per Round';
legend('LEACH','New Algorithm');
xlabel 'Rounds';
ylabel 'Cluster Heads';
hold off
legend('LEACH','New Algorithm');

%%% Plotting Simulation Results "Packet Transfer per Round" 
figure(5)
plot(1:rnd,sim_packet(1:rnd),'-r','Linewidth',2);
hold on
plot(1:rnd,sim_new_packet(1:rnd),'-b','Linewidth',2);
xlim([0 max_rnd]);
title 'Packet Transfer per Round';
legend('LEACH','New Algorithm');
xlabel 'Rounds';
ylabel 'No of Packets';
hold off
legend('LEACH','New Algorithm');

%%% Plotting Simulation Results "Total Energy" 
figure(6)
plot(1:rnd,sim_total_energy(1:rnd),'-r','Linewidth',2);
hold on
plot(1:rnd,sim_new_total_energy(1:rnd),'-b','Linewidth',2);
xlim([0 max_rnd]);
title 'Total Energy Transmitted';
legend('LEACH','New Algorithm');
xlabel 'Rounds';
ylabel 'Energy';
hold off
legend('LEACH','New Algorithm');

%%% Plotting Simulation Results "Average Energy" 
figure(7)
plot(1:rnd,sim_average_energy(1:rnd),'-r','Linewidth',2);
hold on
plot(1:rnd,sim_new_average_energy(1:rnd),'-b','Linewidth',2);
xlim([0 max_rnd]);
title 'Average Energy Transmitted';
legend('LEACH','New Algorithm');
xlabel 'Rounds';
ylabel 'Energy';
hold off
legend('LEACH','New Algorithm');

%%% Plotting Simulation Results "Average Residual Energy" 
figure(8)
plot(1:rnd,sim_average_res_energy(1:rnd),'-r','Linewidth',2);
hold on
plot(1:rnd,sim_new_average_res_energy(1:rnd),'-b','Linewidth',2);
xlim([0 max_rnd]);
ylim([0 inf]);
title 'Average Residual Emergy';
legend('LEACH','New Algorithm');
xlabel 'Rounds';
ylabel 'Energy';
hold off
legend('LEACH','New Algorithm');
