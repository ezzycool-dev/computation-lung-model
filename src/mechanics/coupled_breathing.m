function flow = coupled_breathing(p, lung, L_total, dP_tree_fn, dP_extra_fn)
% COUPLED_BREATHING  Solve the airway-lung flow for one periodic breath.
% Units: V [L], Q_V [L/s], pressures [mmHg], L_total [mmHg s^2/L].
% Two ODEs solved together with ode15s (stiff: airway L/R ~ ms, breath 4 s):
%   dV_A/dt  = Q_V                                   conservation of mass
%   L dQ_V/dt = (P_atm - P_alv) - dP_extra(Q_V) - dP_tree(Q_V)
%                                                    momentum (Newton's 2nd law)
%   P_alv = lung.P_alv_fn(t, V_A, P_pl(t))           downstream lung (linear compliance)
% VT is FIXED: the size of the pleural swing dP_pl is the
% unknown, rescaled until the solved VT equals p.tidal_volume (2-3 passes).
dP_pl = -2.5 * p.cm2mm;          % mmHg, starting guess (Guyton: -5 -> -7.5 cmH2O)
for pass = 1:10
    [flow, VT] = solve_once(p, lung, L_total, dP_tree_fn, dP_extra_fn, dP_pl);
    if abs(VT - p.tidal_volume) < p.VT_tol, break; end
    dP_pl = dP_pl * p.tidal_volume / VT;
end
flow.dP_pl   = dP_pl;           % mmHg, calculated pleural swing
flow.VT      = VT;              % L
flow.passes  = pass;
end


function [flow, VT] = solve_once(p, lung, L_total, dP_tree_fn, dP_extra_fn, dP_pl)
P_pl_fn = @(t) pleural_pressure(t, p, lung, dP_pl);
rhs = @(t, y) [y(2); ...
    (p.p_atm - lung.P_alv_fn(t, y(1), P_pl_fn(t)) ...
     - dP_extra_fn(y(2)) - dP_tree_fn(y(2))) / L_total];

n  = p.n_time;
t_all = linspace(0, p.n_breaths * p.t_cycle, p.n_breaths * (n - 1) + 1)';
opts  = odeset('RelTol', 1e-8, 'AbsTol', [1e-10, 1e-10]);
[t_all, Y] = ode15s(rhs, t_all, [lung.V_FRC; 0], opts);

last = numel(t_all) - n + 1 : numel(t_all);           % final breath
prev = last - (n - 1);                                 % breath before it
flow.periodic_err = max(abs(Y(last, 1) - Y(prev, 1)));   % L

V_A = Y(last, 1);  Q_V = Y(last, 2);
P_pl = zeros(n, 1); P_alv = zeros(n, 1); dQ = zeros(n, 1);
for j = 1:n
    tj = t_all(last(j));
    P_pl(j)  = P_pl_fn(tj);
    P_alv(j) = lung.P_alv_fn(tj, V_A(j), P_pl(j));
    dy       = rhs(tj, Y(last(j), :)');
    dQ(j)    = dy(2);
end

% Flow lags the drive (time constant R*C), so rotate the periodic breath to
% start where flow turns inward: time = 0 is the start of inspiration.
u     = 1:n-1;
start = find(Q_V(u) < 0 & Q_V([2:n-1, 1]) >= 0, 1) + 1;
if isempty(start) || start > n-1, start = 1; end
order = [start:n-1, 1:start-1, start];

flow.time    = t_all(last) - t_all(last(1));
flow.Q_V     = Q_V(order);
flow.dQ_V_dt = dQ(order);
flow.V_A     = V_A(order);
flow.P_alv   = P_alv(order);
flow.P_pl    = P_pl(order);
VT = trapz(flow.time, max(flow.Q_V, 0));
end
