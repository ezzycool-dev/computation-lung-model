function dP = generation_drops(Q_V, p, g)
% GENERATION_DROPS  Resistive pressure drop of each generation in g.
%   Q_V [L/s] (vector, signed). Returns dP [mmHg], numel(Q_V) x numel(g.k).
Q_V = Q_V(:);
dP  = zeros(numel(Q_V), numel(g.k));
for j = 1:numel(Q_V)
    dP(j, :) = resistive_drop(Q_V(j) * 1e-3, g.d, g.l, g.n, g.R_pois_par, ...
                              g.wom_ratio, p.rho_air, p.mu_air, p.use_pedley)' / p.mmHg;
end
end
