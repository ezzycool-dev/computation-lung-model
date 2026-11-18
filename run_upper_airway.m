%% RUN_UPPER_AIRWAY  Solve the Gens 0-4 model, print a summary, save results.
% Works in MATLAB (R2019b or later) and GNU Octave (8+).
%   >> run_upper_airway
root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, 'src')));

p   = upper_airway_params();
out = upper_airway_model(p);
[~, ip] = max(out.Q_V);
b = out.boundaries;

fprintf('\n%s  |  boundary: %s\n', p.version, out.lung_source);
fprintf('------------------------------------------------------------------\n');
fprintf('Breathing     VT %.1f mL, %.0f breaths/min, peak flow %.0f mL/s\n', ...
        out.VT*1e3, 60/p.t_cycle, max(out.Q_V)*1e3);
fprintf('Pressures     alveolar %+.2f mmHg at peak inspiration, pleural %.2f to %.2f mmHg\n', ...
        min(out.P_alv) - p.p_atm, max(out.P_pl), min(out.P_pl));
fprintf('Resistance    nose-larynx %.2f | Gens 0-4 %.3f | Gens 5-16 %.3f mmHg/(L/s)\n', ...
        out.K1, out.R_upper_peak, out.R_low);
fprintf('Drop at peak  nose-larynx %.3f | Gens 0-4 %.4f mmHg\n', ...
        p.p_atm - out.P_trach(ip), out.P_upper(ip,1) - out.P_upper(ip,end));
fprintf('Dead space    %.0f nose-larynx + %.1f Gens 0-4 + %.1f Gens 5-16 = %.0f mL\n', ...
        out.V_extra*1e3, out.V_D_up*1e3, out.V_D_low*1e3, sum(out.V_segs)*1e3);
fprintf('Bohr (CO2)    %.1f mL at the mouth, %.1f mL at the trachea\n', ...
        out.VD_bohr_mouth*1e3, out.VD_bohr_trachea*1e3);
fprintf('Fresh air     reaches the Gen 4 exit %.2f s into inspiration\n', out.t_fresh_gen4);
fprintf('O2 net        %.1f mL/breath at the trachea, %.1f mL at Gen 4/5\n', ...
        b.trachea.mL.net(1), b.gen4_5.mL.net(1));
fprintf('              -> %.0f mL/min O2 uptake, %.0f mL/min CO2 output\n', ...
        b.alveoli.mL.net(1)*60/p.t_cycle, -b.alveoli.mL.net(2)*60/p.t_cycle);
fprintf('------------------------------------------------------------------\n');

export_results(out, fullfile(root, 'results'));
fprintf('Results written to results/  (run scripts/make_figures.py for the figures)\n');
