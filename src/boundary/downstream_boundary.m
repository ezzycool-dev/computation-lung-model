function bnd = downstream_boundary(p)
% DOWNSTREAM_BOUNDARY  What the upper airway sees below Generation 4.
%   bnd = downstream_boundary(p)
%
% This repository models Gens 0-4 in detail. Everything below the Gen 4 exit
% is reduced to the simplest boundary that still lets the airway breathe:
%   * Gens 5-16: plain Hagen-Poiseuille resistance of Weibel's tubes
%     (generations in series, branches in parallel), no gas inertia
%   * Gens 5-16 dead space: carried in the same transport chain, so expired
%     gas reaches Gen 4 with a realistic delay
%   * Alveoli: one linear compliance with fixed alveolar gas (Guyton values)
% The full team model replaces this file with detailed lower-airway,
% respiratory-mechanics and gas-exchange modules.
%
% Returned fields
%   R_total   [mmHg/(L/s)]  lumped Gens 5-16 resistance
%   dP_fn     @(Q_V) -> mmHg resistive drop across Gens 5-16
%   L         [mmHg s^2/L]  inertance (0)
%   vol       [L]           Weibel volume of Gens 5-16 (sets the dead-space split)
%   lung      struct used by coupled_breathing:
%             C_L, V_FRC, P_pl_FRC, shape, P_alv_fn(t, V_A, P_pl), C_alv_fn(t), source

k = p.gens_lower(:);
n = 2.^k;
d = p.weibel_d_cm(k + 1)' / 100;            % m
l = p.weibel_l_cm(k + 1)' / 100;            % m
r = d / 2;
R_si = 8 * p.mu_air * l ./ (pi * r.^4 .* n); % Pa s/m^3 per generation
bnd.R_gen   = R_si * 1e-3 / p.mmHg;          % mmHg/(L/s)
bnd.R_total = sum(bnd.R_gen);
bnd.dP_fn   = @(Q_V) bnd.R_total * Q_V;
bnd.L       = 0;
bnd.vol     = sum(n .* pi .* r.^2 .* l) * 1e3;   % m^3 -> L

% Linear lung: transpulmonary pressure rises by (V_A - V_FRC)/C_L, so the
% alveolar pressure equals atmospheric at FRC.
lung.C_L      = p.C_L;
lung.V_FRC    = p.V_FRC;
lung.P_pl_FRC = p.P_pl_FRC;
lung.shape    = [];                          % half-cosine muscle drive
lung.P_alv_fn = @(t, V_A, P_pl) p.p_atm + P_pl - p.P_pl_FRC + (V_A - p.V_FRC) / p.C_L;
c_alv         = p.pp_alv * p.mmHg / (p.r_gas * p.t_body);
lung.C_alv_fn = @(t) c_alv;
lung.source   = 'linear-compliance lung (Guyton values)';
bnd.lung = lung;
end
