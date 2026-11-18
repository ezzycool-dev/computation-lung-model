function [R_pois_par, wom_ratio, L_par] = womersley_terms(d, l, n, rho, mu, omega)
% WOMERSLEY_TERMS  Flow-independent impedance of each generation (SI).
%   d, l [m], n tubes in parallel; rho [kg/m^3], mu [Pa s], omega [rad/s].
% Womersley's oscillatory solution of the Navier-Stokes equations in a rigid
% tube (Bessel functions): Z = i w rho L / (pi r^2 (1 - F10)),
%   F10 = 2 J1(Lam) / (Lam J0(Lam)),  Lam = i^(3/2) Wo,  Wo = r sqrt(rho w / mu).
% Outputs (per generation, n tubes in parallel):
%   R_pois_par [Pa s/m^3]  Hagen-Poiseuille resistance 8 mu L / (pi r^4) / n
%   wom_ratio  [-]         Re(Z) / R_Poiseuille
%   L_par      [Pa s^2/m^3] inertance Im(Z)/omega / n
r          = d / 2;
Wo         = r .* sqrt(rho * omega / mu);
R_pois     = 8 * mu .* l ./ (pi * r.^4);
Lam        = Wo * 1i^(3/2);
F10        = 2 * besselj(1, Lam) ./ (Lam .* besselj(0, Lam));
Z          = 1i * omega * rho .* l ./ (pi * r.^2) ./ (1 - F10);
wom_ratio  = real(Z) ./ R_pois;
R_pois_par = R_pois ./ n;
L_par      = imag(Z) / omega ./ n;
end
