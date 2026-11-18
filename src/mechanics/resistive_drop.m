function dP = resistive_drop(Q, d, l, n, R_pois_par, wom_ratio, rho, mu, use_pedley)
% RESISTIVE_DROP  Viscous pressure drop across each generation (SI).
%   Q [m^3/s, signed], geometry in m. Returns dP [Pa] per generation, sign of Q.
% R = R_Poiseuille * max(Womersley ratio, Pedley ratio, 1)
%   Pedley, Schroter & Sudlow (1970): R/R_Pois = 0.327 sqrt(Re d / L)
%   Re = rho u d / mu, u = |Q| / total area of the generation (continuity)
A   = n .* pi .* (d/2).^2;
Re  = rho * (abs(Q) ./ A) .* d / mu;
ped = 0.327 * sqrt(Re .* d ./ l);
if use_pedley
    ratio = max(max(wom_ratio, ped), 1);
else
    ratio = wom_ratio;
end
dP = R_pois_par .* ratio * Q;
end
