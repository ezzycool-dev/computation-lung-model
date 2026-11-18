function P = upper_profile(time, Q_V, dQ_V_dt, P_top, p, g)
% UPPER_PROFILE  Node pressures along Gens 0-4 over the breath.
%   time    [s]      N x 1
%   Q_V     [L/s]    N x 1
%   dQ_V_dt [L/s^2]  N x 1
%   P_top   [mmHg]   N x 1 absolute pressure at the trachea entry
%   P       [mmHg]   N x (5+1) absolute: column 1 = Gen 0 entry, end = Gen 4 exit
% Momentum balance per generation (Newton's 2nd law on the gas column):
%   P_k - P_{k+1} = R_k(Q) Q + L_k dQ/dt
drops = generation_drops(Q_V, p, g) + dQ_V_dt(:) * g.L_gen';
P = P_top(:) - [zeros(numel(time), 1), cumsum(drops, 2)];
end
