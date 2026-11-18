function [K1, cal] = calibrate_extrathoracic(p, dP_tree_fn)
% CALIBRATE_EXTRATHORACIC  Nose-to-larynx laminar coefficient K1 [mmHg/(L/s)].
% From Guyton Ch. 37: ~1 cmH2O alveolar pressure "is enough to pull 0.5 liter
% of air into the lungs in the 2 seconds required for normal quiet inspiration".
% 1) Peak flow of a half-sine inspiration moving VT_ref in T_insp_ref:
%       Q_pk = pi * VT_ref / (2 * T_insp_ref)
% 2) Whole-airway resistance at peak flow (Ohm's-law analogue, Sys Phys 6.2.3):
%       R_total = dP_alv_ref / Q_pk
% 3) Resistances in series add (Rohrer, Sys Phys Eqn 6.2.10), so the part not
%    explained by Gens 0-16 belongs to the nose/mouth/pharynx/larynx:
%       K1 = R_total - R_tree(Q_pk)
% dP_tree_fn(Q_V) must return the Gens 0-16 resistive drop [mmHg] for Q_V [L/s].
cal.Q_pk    = pi * p.VT_ref / (2 * p.T_insp_ref);       % L/s
cal.R_total = p.dP_alv_ref / cal.Q_pk;                 % mmHg/(L/s)
cal.R_tree  = dP_tree_fn(cal.Q_pk) / cal.Q_pk;         % mmHg/(L/s)
K1 = cal.R_total - cal.R_tree;
if K1 < 0
    warning('Airway tree resistance exceeds the Guyton total; K1 set to 0.');
    K1 = 0;
end
cal.K1 = K1;
end
