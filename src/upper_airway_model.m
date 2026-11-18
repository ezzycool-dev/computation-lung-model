function out = upper_airway_model(p, bnd)
% UPPER_AIRWAY_MODEL  Airflow, pressure and dead-space gas transport in
% airway Generations 0-4 (trachea to segmental bronchi).
%   out = upper_airway_model()            % default parameters
%   out = upper_airway_model(p)           % p = upper_airway_params(), edited
%   out = upper_airway_model(p, bnd)      % bnd = downstream_boundary(p), edited
%
% Solve order
%   1  Weibel geometry and impedance of Gens 0-4
%   2  airway-tree resistance (Gens 0-4 detailed + Gens 5-16 boundary)
%   3  nose-to-larynx resistance K1, calibrated to Guyton's quiet breath
%   4  coupled airway-lung flow with the tidal volume held at 500 mL
%   5  pressure at every generation boundary
%   6  gas transport through a chain of well-mixed segments
%   7  gas crossing each boundary, Bohr dead space, packaging
% Units: mmHg, L/s, L, s, mol/m^3 (see upper_airway_params).

if nargin < 1 || isempty(p),   p   = upper_airway_params(); end
if nargin < 2 || isempty(bnd), bnd = downstream_boundary(p); end

%% 1-2. Geometry and resistance
gu = upper_geometry(p);
dP_tree_fn = @(Q) upper_dP(Q, p, gu) + bnd.dP_fn(Q);            % mmHg

%% 3. Nose-to-larynx resistance (Rohrer form, linear term calibrated)
[K1, cal]   = calibrate_extrathoracic(p, dP_tree_fn);
dP_extra_fn = @(Q) K1 * Q + p.K2 * Q .* abs(Q);                % mmHg

%% 4. Coupled breathing: one flow for the whole airway
flow = coupled_breathing(p, bnd.lung, gu.L + bnd.L, dP_tree_fn, dP_extra_fn);
t = flow.time;  Q = flow.Q_V;

%% 5. Pressures: mouth -> trachea -> each Gen 0-4 boundary -> Gen 4 exit
P_trach = p.p_atm - dP_extra_fn(Q);
P_upper = upper_profile(t, Q, flow.dQ_V_dt, P_trach, p, gu);   % N x 6
P_gen4  = P_upper(:, end);

%% 6. Dead-space split and gas transport
% The 150 mL dead space = 72 mL nose-larynx + intrathoracic part. The
% intrathoracic part is shared between Gens 0-4 and Gens 5-16 in proportion
% to Weibel's airway volumes.
V_intra = p.dead_space_vol - p.V_extra;
frac_up = sum(gu.vol) * 1e3 / (sum(gu.vol) * 1e3 + bnd.vol);
V_D_up  = V_intra * frac_up;                                   % L, Gens 0-4
V_D_low = V_intra - V_D_up;                                    % L, Gens 5-16
V_seg   = V_D_up / p.n_seg_upper;
n_low   = max(1, round(V_D_low / V_seg));
if p.track_extrathoracic
    n_extra = max(1, round(p.V_extra / V_seg));
else
    n_extra = 0;
end
V_segs = [repmat(p.V_extra / max(n_extra, 1), n_extra, 1);
          repmat(V_seg, p.n_seg_upper, 1);
          repmat(V_D_low / n_low, n_low, 1)];
i_tr = n_extra + 1;                 % first segment of the trachea (Gen 0)
i_4  = n_extra + p.n_seg_upper;     % last segment of Gen 4

c_init = repmat(bnd.lung.C_alv_fn(0), numel(V_segs), 1);
gas = segment_transport(t, Q, V_segs, p.c_trach, bnd.lung.C_alv_fn, p.n_breaths, c_init);
C = gas.C;                          % N x n_seg x 4

%% 7. Gas crossing each boundary (upwind: the segment the gas leaves)
ins  = Q >= 0;
segC = @(j) reshape(C(:, j, :), numel(t), 4);
C_alv_t = cell2mat(arrayfun(@(tt) bnd.lung.C_alv_fn(tt), t, 'UniformOutput', false));
bounds = struct();
bounds.mouth   = face(repmat(p.c_trach, numel(t), 1), segC(1), ins);
if n_extra > 0
    bounds.trachea = face(segC(i_tr - 1), segC(i_tr), ins);
else
    bounds.trachea = bounds.mouth;
end
bounds.gen4_5  = face(segC(i_4), segC(i_4 + 1), ins);
bounds.alveoli = face(segC(numel(V_segs)), C_alv_t, ins);
names = fieldnames(bounds);
for b = 1:numel(names)
    bounds.(names{b}).mL = crossing_volumes(t, Q, bounds.(names{b}).C);
end

%% Bohr dead space from the CO2 trace at a boundary (capnography)
VD_mouth   = bohr_dead_space(t, Q, bounds.mouth.C(:, 2), p, bnd);
VD_trachea = bohr_dead_space(t, Q, bounds.trachea.C(:, 2), p, bnd);

% Inspired volume at which fresh air first reaches the Gen 4 exit
cumV = cumtrapz(t, max(Q, 0));
i_in = 1:find(Q(2:end) < 0, 1);
t_fresh_gen4 = interp1(cumV(i_in), t(i_in), V_D_up + n_extra * (p.V_extra / max(n_extra, 1)));

%% Per-generation table at peak inspiratory flow
[~, ip] = max(Q);
gen = generation_table(Q(ip), p, gu);

%% Package
out.params = p;
out.time = t;                 out.Q_V = Q;              out.dQ_V_dt = flow.dQ_V_dt;
out.VT = flow.VT;             out.V_A = flow.V_A;
out.P_alv = flow.P_alv;       out.P_pl = flow.P_pl;     out.dP_pl = flow.dP_pl;
out.P_mouth = p.p_atm * ones(size(t));
out.P_trach = P_trach;        out.P_upper = P_upper;    out.P_gen4 = P_gen4;
out.K1 = K1;                  out.calibration = cal;    out.R_low = bnd.R_total;
out.R_upper_peak = sum(gen.R_mmHg_per_Lps);
out.geometry = gu;            out.gen_table = gen;      out.i_peak = ip;
out.V_extra = p.V_extra;      out.V_D_up = V_D_up;      out.V_D_low = V_D_low;
out.V_segs = V_segs;          out.i_trachea = i_tr;     out.i_gen4 = i_4;
out.n_extra = n_extra;
out.C_chain = C;              out.boundaries = bounds;
out.VD_bohr_mouth = VD_mouth; out.VD_bohr_trachea = VD_trachea;
out.t_fresh_gen4 = t_fresh_gen4;
out.flow_periodic_err = flow.periodic_err;
out.gas_periodic_err  = max(abs(gas.y_end_last - gas.y_start_last));
out.lung_source = bnd.lung.source;
end


function f = face(C_above, C_below, ins)
% Gas crossing a face: from the segment above on inspiration, from the one below on expiration.
f.C = C_below;
f.C(ins, :) = C_above(ins, :);
end


function v = crossing_volumes(t, Q, Cf)
% Partial volume (mL) of each gas crossing a face, inspired and expired.
frac = Cf ./ sum(Cf, 2);
qin  = max(Q, 0);   qout = max(-Q, 0);
v.inspired = trapz(t, qin  .* frac) * 1e3;      % 1 x 4, mL
v.expired  = trapz(t, qout .* frac) * 1e3;      % 1 x 4, mL
v.net      = v.inspired - v.expired;            % 1 x 4, mL (+ = toward the alveoli)
end


function VD = bohr_dead_space(t, Q, c_co2, p, bnd)
% Bohr: VD/VT = (P_A,CO2 - P_E,CO2)/(P_A,CO2 - P_I,CO2), with P_E the
% flow-weighted mean expired CO2 at the measuring point.
to_mm = @(c) c * p.r_gas * p.t_body / p.mmHg;
q  = max(-Q, 0);
PE = trapz(t, q .* to_mm(c_co2)) / trapz(t, q);
ca = bnd.lung.C_alv_fn(0);
PA = to_mm(ca(2));
PI = p.pp_trach(2);
VD = trapz(t, q) * (PA - PE) / (PA - PI);
end


function g = generation_table(Qpk, p, gu)
% Reynolds, Womersley and Pedley numbers, resistance and pressure drop of
% each generation at peak inspiratory flow.
Qs   = Qpk * 1e-3;                                  % m^3/s
u    = Qs ./ gu.area;                               % m/s, mean velocity
Re   = p.rho_air * u .* gu.d / p.mu_air;
ped  = 0.327 * sqrt(Re .* gu.d ./ gu.l);
ratio = max(max(gu.wom_ratio, ped), 1);
if ~p.use_pedley, ratio = gu.wom_ratio; end
R    = gu.R_pois_par .* ratio * 1e-3 / p.mmHg;      % mmHg/(L/s)
g.gen = gu.k;           g.n = gu.n;
g.d_cm = gu.d * 100;    g.l_cm = gu.l * 100;
g.area_cm2 = gu.area * 1e4;   g.vol_mL = gu.vol * 1e6;
g.velocity_cm_s = u * 100;
g.Re = Re;              g.Wo = gu.Wo;
g.womersley_ratio = gu.wom_ratio;   g.pedley_ratio = ped;   g.R_ratio = ratio;
g.R_poiseuille = gu.R_pois_par * 1e-3 / p.mmHg;
g.R_mmHg_per_Lps = R;
g.dP_peak_mmHg = R * Qpk;
g.share_pct = 100 * R / sum(R);
end
