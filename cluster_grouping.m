function [SN] = cluster_grouping(SN,CL, CLheads)
%CLUSTER_GROUPING Grouping Nodes into Clusters
%   This function group nodes in the wireless sensor network (WSN) into the
%   clusters based on their distance from the elected cluster heads
%
%   INPUT PARAMETERS
%   SN - all sensors nodes (including routing routes)
%   CL - all cluster nodes
%   CLheads - number of cluster heads elected.
%
%   OUTPUT PARAMETERS
%   SN - all sensors nodes (including routing routes)

for i=1:length(SN.n)
    if  (SN.n(i).role == 'N') && (SN.n(i).E > 0) && CLheads > 0 % if node is normal

        dist = zeros(1, length(CL.n));

        for j = 1:length(CL.n)
            % distance between the sensor node and each cluster head
            dist(j)=sqrt( (CL.n(j).x-SN.n(i).x)^2 + (CL.n(j).y-SN.n(i).y)^2 );
        end

        [~, index] = min( dist(:) ); % finds the minimum distance of node to CH
        SN.n(i).cluster = index; % assigns node to the cluster
        SN.n(i).dnc = dist(index); % assigns the distance of node to CH
        SN.n(i).chid = CL.n(index).id;
    end
end

end

