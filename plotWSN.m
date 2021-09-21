function plotWSN(SN, dims, sim_name)
%PLOTWSN Graphical Representation of Wireless Sensor Network
%   This function provides of the graphical representation of a Wireless 
%   Sensor Network (WSN). 
%
%   INPUT PARAMETERS
%   SN - all sensors nodes (including routing nodes)
%   dims - container of the dimensions of the WSN plot extremes and the
%           base station point. outputs: x_min, x_min, y_min, x_max, y_max, 
%           bs_x, bs_y
%   rn_circle - a boolean for plotting routing node line. DEFAULT: false
%   sim_name - the name of the wireless network simulation. 
%               Default: 'LEACH'.
%
%   OUTPUT: the graphical representation of the WSN

if nargin < 3
    sim_name = 'LEACH';
end

%{
if (~ischar(sim_name)) && (~isstring(sim_name))
    error('The details passed for ''the sim_name'' argument is not a string. \nString is required')
end
%}

hold on

plot( dims('x_min'),dims('y_min'),dims('x_max'),dims('y_max') )
for i=1:length(SN.n)
    if strcmp( SN.n(i).role, 'N' )
        plot(SN.n(i).x,SN.n(i).y,'om');
    elseif strcmp( SN.n(i).role, 'R' )
        plot(SN.n(i).x,SN.n(i).y, 'or','MarkerFaceColor','b','MarkerEdgeColor','r','MarkerSize',8);
    end
end

plot( dims('bs_x'), dims('bs_y'), '*k','LineWidth',3,'MarkerSize',10);
title ({sim_name; 'Wireless Sensor Network';})
xlabel '(m)';
ylabel '(m)';

if dims('rn_dist') ~= 0
    theta = linspace(0,2*pi);
    xcir = dims('rn_dist')*cos(theta) + dims('bs_x');
    ycir = dims('rn_dist')*sin(theta) + dims('bs_y'); 
    
    plot(xcir,ycir,'-b','LineWidth',3);
end

end

