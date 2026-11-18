function export_results(out, folder)
% EXPORT_RESULTS  Write the model's results as plain CSV files (read by
% scripts/make_figures.py and easy to open in Excel).
%   export_results(out, 'results')
if nargin < 2, folder = 'results'; end
if ~exist(folder, 'dir'), mkdir(folder); end
p  = out.params;
mm = @(c) c * p.r_gas * p.t_body / p.mmHg;          % mol/m^3 -> mmHg

%% 1. Time series over one breath
b = out.boundaries;
ts = [out.time, out.Q_V, out.V_A - out.V_A(1), out.P_alv - p.p_atm, out.P_pl, ...
      p.p_atm - out.P_trach, p.p_atm - out.P_upper, ...
      mm(b.mouth.C(:, 2)), mm(b.trachea.C(:, 2)), mm(b.gen4_5.C(:, 2)), mm(b.alveoli.C(:, 2)), ...
      mm(b.mouth.C(:, 1)), mm(b.trachea.C(:, 1)), mm(b.gen4_5.C(:, 1)), mm(b.alveoli.C(:, 1))];
hdr = {'t_s', 'Q_Lps', 'dV_L', 'P_alv_gauge_mmHg', 'P_pl_gauge_mmHg', ...
       'drop_trachea_mmHg', 'drop_gen0_in_mmHg', 'drop_gen1_in_mmHg', 'drop_gen2_in_mmHg', ...
       'drop_gen3_in_mmHg', 'drop_gen4_in_mmHg', 'drop_gen4_out_mmHg', ...
       'PCO2_mouth', 'PCO2_trachea', 'PCO2_gen4_5', 'PCO2_alveoli', ...
       'PO2_mouth', 'PO2_trachea', 'PO2_gen4_5', 'PO2_alveoli'};
write_csv(fullfile(folder, 'timeseries.csv'), hdr, ts);

%% 2. Space-time map of CO2 along the chain (mouth -> Gen 4/5 handoff)
idx  = 1:(out.i_gen4 + 1);
vcum = cumsum(out.V_segs) - out.V_segs / 2;          % L, centre of each segment
PCO2 = mm(out.C_chain(:, idx, 2));
write_csv(fullfile(folder, 'co2_map.csv'), ...
          [{'t_s'}, arrayfun(@(v) sprintf('V%.2fmL', v), vcum(idx)' * 1e3, 'UniformOutput', false)], ...
          [out.time, PCO2]);
write_csv(fullfile(folder, 'co2_map_axis.csv'), {'segment', 'volume_from_mouth_mL', 'region'}, ...
          [idx', vcum(idx) * 1e3, region_code(idx, out)']);

%% 3. Generation table at peak flow
g = out.gen_table;
write_csv(fullfile(folder, 'generations.csv'), ...
    {'gen', 'n_tubes', 'd_cm', 'l_cm', 'area_cm2', 'vol_mL', 'velocity_cm_s', 'Re', 'Wo', ...
     'womersley_ratio', 'pedley_ratio', 'R_ratio', 'R_poiseuille', 'R_mmHg_per_Lps', ...
     'dP_peak_mmHg', 'share_pct'}, ...
    [g.gen, g.n, g.d_cm, g.l_cm, g.area_cm2, g.vol_mL, g.velocity_cm_s, g.Re, g.Wo, ...
     g.womersley_ratio, g.pedley_ratio, g.R_ratio, g.R_poiseuille, g.R_mmHg_per_Lps, ...
     g.dP_peak_mmHg, g.share_pct]);

%% 4. Gas crossing each boundary per breath (mL)
names = {'mouth', 'trachea', 'gen4_5', 'alveoli'};
M = zeros(numel(names), 12);
for k = 1:numel(names)
    v = b.(names{k}).mL;
    M(k, :) = [v.inspired, v.expired, v.net];
end
fid = fopen(fullfile(folder, 'boundary_gas.csv'), 'w');
fprintf(fid, 'boundary,O2_in,CO2_in,N2_in,H2O_in,O2_out,CO2_out,N2_out,H2O_out,O2_net,CO2_net,N2_net,H2O_net\n');
for k = 1:numel(names)
    fprintf(fid, '%s', names{k}); fprintf(fid, ',%.4f', M(k, :)); fprintf(fid, '\n');
end
fclose(fid);

%% 5. Headline numbers
[~, ip] = max(out.Q_V);
s = {
  'tidal_volume_mL',            out.VT * 1e3
  'peak_flow_mL_s',             max(out.Q_V) * 1e3
  'breaths_per_min',            60 / p.t_cycle
  'K1_nose_larynx_mmHg_per_Lps', out.K1
  'R_gen0_4_peak_mmHg_per_Lps', out.R_upper_peak
  'R_gen5_16_mmHg_per_Lps',     out.R_low
  'drop_nose_larynx_peak_mmHg', p.p_atm - out.P_trach(ip)
  'drop_gen0_4_peak_mmHg',      out.P_upper(ip, 1) - out.P_upper(ip, end)
  'P_alv_min_gauge_mmHg',       min(out.P_alv) - p.p_atm
  'P_pl_min_mmHg',              min(out.P_pl)
  'P_pl_max_mmHg',              max(out.P_pl)
  'dead_space_nose_larynx_mL',  out.V_extra * 1e3
  'dead_space_gen0_4_mL',       out.V_D_up * 1e3
  'dead_space_gen5_16_mL',      out.V_D_low * 1e3
  'segments_gen0_4',            p.n_seg_upper
  'segment_volume_mL',          out.V_D_up / p.n_seg_upper * 1e3
  'bohr_dead_space_mouth_mL',   out.VD_bohr_mouth * 1e3
  'bohr_dead_space_trachea_mL', out.VD_bohr_trachea * 1e3
  't_fresh_air_reaches_gen4_s', out.t_fresh_gen4
  'O2_uptake_mL_min',           b.alveoli.mL.net(1) * 60 / p.t_cycle
  'CO2_output_mL_min',          -b.alveoli.mL.net(2) * 60 / p.t_cycle
  'alveolar_ventilation_L_min', (out.VT - p.dead_space_vol) * 60 / p.t_cycle
  'flow_periodic_err_L',        out.flow_periodic_err
  'gas_periodic_err_mol_m3',    out.gas_periodic_err};
fid = fopen(fullfile(folder, 'summary.csv'), 'w');
fprintf(fid, 'quantity,value\n');
for k = 1:size(s, 1)
    fprintf(fid, '%s,%.6g\n', s{k, 1}, s{k, 2});
end
fclose(fid);
end


function code = region_code(idx, out)
% 0 = nose-larynx, 1 = Gens 0-4, 2 = Gens 5-16
code = 2 * ones(size(idx));
code(idx < out.i_trachea) = 0;
code(idx >= out.i_trachea & idx <= out.i_gen4) = 1;
end


function write_csv(fname, header, M)
fid = fopen(fname, 'w');
fprintf(fid, '%s', header{1});
fprintf(fid, ',%s', header{2:end});
fprintf(fid, '\n');
fmt = [repmat('%.8g,', 1, size(M, 2) - 1), '%.8g\n'];
fprintf(fid, fmt, M');
fclose(fid);
end
