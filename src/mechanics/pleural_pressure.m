function P = pleural_pressure(t, p, lung, dP_pl)
% PLEURAL_PRESSURE  Muscle drive [mmHg gauge] = P_pl_FRC + dP_pl * s(t).
% s(t) goes from 0 (FRC) to 1 (end of inspiration): lung.shape if
% supplied, otherwise a half-cosine rise over t_in and fall over t_out.
% dP_pl (the size of the swing) is calculated by coupled_breathing so that
% the solved tidal volume equals p.tidal_volume.
if isempty(lung.shape)
    tc = mod(t, p.t_cycle);
    if tc <= p.t_in
        s = (1 - cos(pi * tc / p.t_in)) / 2;
    else
        s = (1 + cos(pi * (tc - p.t_in) / p.t_out)) / 2;
    end
else
    tp = mod(t, lung.shape(end, 1));
    s  = interp1(lung.shape(:, 1), lung.shape(:, 2), tp, 'linear');
end
P = lung.P_pl_FRC + dP_pl * s;
end
