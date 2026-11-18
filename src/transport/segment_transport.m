function gas = segment_transport(time, Q_V, V_segs, c_top, C_bottom_fn, n_breaths, c_init)
% SEGMENT_TRANSPORT  Gas along one chain of well-mixed segments.
%   time, Q_V   shared time grid [s] and flow [L/s] (+ = toward the alveoli)
%   V_segs      segment volumes [L], segment 1 = trachea entry, end = Gen-16 exit
%   c_top       1 x 4 gas entering at the trachea on inspiration [mol/m^3]
%   C_bottom_fn @(t) -> 1 x 4 gas entering from the alveoli on expiration
%   c_init      n_seg x 4 initial contents (end-expiration)
% Species mass balance on each well-mixed segment (in - out, tanks in series):
%   inspiration: V_j dC_j/dt = Q_V (C_{j-1} - C_j),  C_0 = c_top
%   expiration:  V_j dC_j/dt = |Q_V| (C_{j+1} - C_j), C_{n+1} = C_bottom(t)
% As segments shrink this is the advection equation dC/dt + Q dC/dV = 0.
% Solved with ode45; repeated n_breaths times to reach periodic steady state.
% Output gas.C is numel(time) x n_seg x 4 for the final breath.
V_segs = V_segs(:);
n_seg  = numel(V_segs);
n_sp   = size(c_init, 2);
opts   = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'MaxStep', min(diff(time)));
rhs    = @(t, y) segment_rhs(t, y, time, Q_V, V_segs, n_sp, c_top, C_bottom_fn);

y0 = c_init(:);
gas.y_start_prev = y0;
for b = 1:n_breaths
    gas.y_start_prev = y0;
    [~, Y] = ode45(rhs, time, y0, opts);
    y0 = Y(end, :)';
end
gas.y_start_last = gas.y_start_prev;
gas.y_end_last   = Y(end, :)';
gas.C = reshape(Y, numel(time), n_seg, n_sp);
end


function dydt = segment_rhs(t, y, time, Q_V, V_segs, n_sp, c_top, C_bottom_fn)
q = interp1(time, Q_V, t, 'linear', 'extrap');
C = reshape(y, numel(V_segs), n_sp);
if q >= 0
    upstream = [c_top; C(1:end-1, :)];
else
    upstream = [C(2:end, :); C_bottom_fn(t)];
end
dCdt = abs(q) ./ V_segs .* (upstream - C);
dydt = dCdt(:);
end
