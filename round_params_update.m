function [round_params, stability_period_check, lifetime_check] = round_params_update(SN, CLheads, round_params, rn_ids, round, rounds, stability_period_check, lifetime_check)
%ROUND_PARAMS_UPDATE Update the Simulation Parameters during a round
%   This function is used to update all the parameters used in  gathering
%   data for the analytics of the wireless network sensor (WSN).
%
%   INPUT PARAMETERS
%   SN - all sensors nodes (including routing routes)
%   CLheads - number of cluster heads elected.
%   round_params - container of the parameters used to measure the
%                   performance of the simulation in a round. The params
%                   are: 'dead nodes', 'operating nodes', 'total energy', 
%                   'packets', 'stability period', 'lifetime', 
%                   'stability period round', 'lifetime round'.
%   rn_ids - ids of all sensor nodes used for routing
%   round - the current round in the simulation.
%   rounds - the total number of rounds in the simualtion.
%   stability_period_check - boolean indicating the active search of the
%                               stability period metric.
%   lifetime_check - boolean indication the active search of the lifetime
%                       period metric.
%
%   OUTPUT PARAMETERS
%   round_params - container of the parameters used to measure the
%                   performance of the simulation in a round. The params
%                   are: 'dead nodes', 'operating nodes', 'total energy', 
%                   'packets', 'stability period', 'lifetime', 
%                   'stability period round', 'lifetime round'.
%   stability_period_check - boolean indicating the active search of the
%                               stability period metric.
%   lifetime_check - boolean indication the active search of the lifetime
%                       period metric.

if stability_period_check
    if round_params('operating nodes') < length(SN.n) - length(rn_ids)
        round_params('stability period') = toc;
        round_params('stability period round') = round;
        stability_period_check = false;
    elseif round == rounds
        round_params('stability period') = toc;
        round_params('stability period round') = round;
    end
end

if lifetime_check
    if round_params('operating nodes') == length(rn_ids)
        round_params('lifetime') = toc;
        round_params('lifetime round') = round;
        lifetime_check = false;
    elseif round == rounds
        round_params('lifetime') = toc;
        round_params('lifetime round') = round;
    end
end

round_params('cluster heads') = CLheads;

end

