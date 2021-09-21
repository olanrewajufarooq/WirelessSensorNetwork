%% Network Establishment Parameters 

% Number of simulations
sim_n=3;

% Size of data package %
k=800; % units in bits

% Suggested percentage of cluster head %
p=0.05; % a 5 percent of the total amount of nodes used in the network is proposed to give good results

% Maximum Number of Rounds %
max_rnd = 6000;


%% Parameters for LEACH Algorithm

% Number of Nodes in the field %
n=104;


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

        % Reseting Previous Amount Of Energy Consumed In the Network on the Previous Round %
        energy=0;
        % Reseting Previous Amount Of Energy Consumed In the Network on the Previous Round %
        new_energy=0;
        
        
       %%%%%% Energy Dissipation for normal nodes %%%%%%


        % LEACH %

        

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
