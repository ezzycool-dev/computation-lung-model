function g = weibel_geometry(p, gens)
% WEIBEL_GEOMETRY  Geometry and flow-independent impedance for a set of
% Weibel generations (symmetric dichotomy, n = 2^k tubes per generation).
% Returns struct g:
%   k, n          generation numbers and tube counts
%   d, l          diameter and length of one tube [m]
%   area, vol     total cross-section [m^2] and volume [m^3] per generation
%   R_pois_par, wom_ratio   resistance terms (SI, see womersley_terms)
%   L_gen         inertance per generation [mmHg s^2/L];  L = sum(L_gen)
%   Wo            Womersley number per generation
omega  = 2*pi / p.t_cycle;
g.k    = gens(:);
g.n    = 2.^g.k;
g.d    = p.weibel_d_cm(g.k + 1)' / 100;
g.l    = p.weibel_l_cm(g.k + 1)' / 100;
g.area = g.n .* pi .* (g.d/2).^2;
g.vol  = g.area .* g.l;
[g.R_pois_par, g.wom_ratio, L_par] = womersley_terms(g.d, g.l, g.n, ...
                                         p.rho_air, p.mu_air, omega);
g.L_gen = L_par * 1e-3 / p.mmHg;      % Pa s^2/m^3 -> mmHg s^2/L
g.L     = sum(g.L_gen);
g.Wo    = (g.d/2) .* sqrt(p.rho_air * omega / p.mu_air);
end
