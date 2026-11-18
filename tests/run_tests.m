function results = run_tests()
% RUN_TESTS  Physical and numerical checks for the upper-airway model.
%   results = run_tests();    prints PASS/FAIL per check and returns a struct array.
% Each check compares the model with a textbook value, an exact conservation
% law, or itself at a finer resolution.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')));

p   = upper_airway_params();
bnd = downstream_boundary(p);
o   = upper_airway_model(p, bnd);
b   = o.boundaries;
[~, ip] = max(o.Q_V);
R = {};

%% Mechanics
R(end+1, :) = check('Tidal volume held at 500 mL', ...
    abs(o.VT - p.tidal_volume) < 1e-4, sprintf('%.2f mL', o.VT*1e3));

a = min(o.P_alv) - p.p_atm;
R(end+1, :) = check('Alveolar P at peak inspiration = -1 cmH2O (Guyton)', ...
    abs(a + p.dP_alv_ref) < 0.05, sprintf('%.3f mmHg (target %.3f)', a, -p.dP_alv_ref));

e = max(abs(o.P_gen4 - bnd.dP_fn(o.Q_V) - o.P_alv));
R(end+1, :) = check('Pressure path mouth -> Gen 4 -> alveoli closes', ...
    e < 1 / p.mmHg, sprintf('%.1e mmHg', e));

ratio = o.R_upper_peak / o.R_low;
R(end+1, :) = check('Gens 0-4 carry ~2x the resistance of Gens 5-16', ...
    ratio > 1.5 && ratio < 3, sprintf('%.2fx', ratio));

R(end+1, :) = check('Laminar flow at rest: Re < 2000 in every generation (Feher)', ...
    all(o.gen_table.Re < 2000), sprintf('max Re = %.0f (trachea)', max(o.gen_table.Re)));

%% Gas transport
dO2  = [b.trachea.mL.net(1), b.gen4_5.mL.net(1), b.alveoli.mL.net(1)];
dCO2 = [b.trachea.mL.net(2), b.gen4_5.mL.net(2), b.alveoli.mL.net(2)];
R(end+1, :) = check('Mass conservation: same net O2 and CO2 at every boundary', ...
    range(dO2) < 0.05 && range(dCO2) < 0.05, ...
    sprintf('O2 %.2f/%.2f/%.2f mL, CO2 %.2f/%.2f/%.2f mL', dO2, dCO2));

VA   = o.VT - sum(o.V_segs);
fick = VA * 1e3 * (p.pp_trach(1:2) - p.pp_alv(1:2)) / p.p_atm;   % mL per breath
net  = b.alveoli.mL.net(1:2);
R(end+1, :) = check('Fick: net O2/CO2 = alveolar ventilation x partial-pressure gap', ...
    all(abs(net - fick) ./ abs(fick) < 0.01), ...
    sprintf('O2 %.2f vs %.2f mL, CO2 %.2f vs %.2f mL', net(1), fick(1), net(2), fick(2)));

R(end+1, :) = check('Bohr dead space at the mouth = 150 mL (Guyton)', ...
    abs(o.VD_bohr_mouth - p.dead_space_vol) < 0.002, sprintf('%.1f mL', o.VD_bohr_mouth*1e3));

vi = sum(o.V_segs) - o.V_extra;
R(end+1, :) = check('Bohr dead space at the trachea = intrathoracic volume', ...
    abs(o.VD_bohr_trachea - vi) < 0.002, sprintf('%.1f mL (chain %.1f mL)', o.VD_bohr_trachea*1e3, vi*1e3));

C = reshape(o.C_chain, [], 4);
eg = check_ideal_gas(C, p.p_atm, p);
R(end+1, :) = check('Ideal gas: species sum to P/RT everywhere', ...
    eg < 0.01, sprintf('%.1e mol/m^3', eg));

%% Numerics
R(end+1, :) = check('Periodic steady state (flow and gas)', ...
    o.flow_periodic_err < 1e-4 && o.gas_periodic_err < 1e-3, ...
    sprintf('flow %.1e L, gas %.1e mol/m^3', o.flow_periodic_err, o.gas_periodic_err));

p2 = p;  p2.n_seg_upper = 2 * p.n_seg_upper;
o2 = upper_airway_model(p2, downstream_boundary(p2));
R(end+1, :) = check('Converged: doubling the segment count moves dead space < 1 mL', ...
    abs(o2.VD_bohr_mouth - o.VD_bohr_mouth) < 0.001, ...
    sprintf('%.2f vs %.2f mL', o2.VD_bohr_mouth*1e3, o.VD_bohr_mouth*1e3));

drop = o.P_upper(ip, 1) - o.P_upper(ip, end);
R(end+1, :) = check('Regression vs. earlier version (peak flow, Gens 0-4 drop)', ...
    abs(max(o.Q_V)*1e3 - 385.8)/385.8 < 0.01 && abs(drop - 0.0738)/0.0738 < 0.02, ...
    sprintf('%.1f mL/s, %.4f mmHg', max(o.Q_V)*1e3, drop));

n_pass = sum(cell2mat(R(:, 2)));
fprintf('\n%d of %d checks passed\n', n_pass, size(R, 1));
results = cell2struct(R, {'name', 'pass', 'value'}, 2);
end


function row = check(name, ok, value)
persistent k
if isempty(k), k = 0; end
k = k + 1;
tag = 'FAIL'; if ok, tag = 'PASS'; end
fprintf('%s  %2d  %-66s %s\n', tag, k, name, value);
row = {name, logical(ok), value};
end


function r = range(x)
r = max(x) - min(x);
end
