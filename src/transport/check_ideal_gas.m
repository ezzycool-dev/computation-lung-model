function err = check_ideal_gas(C, P_mmHg, p)
% CHECK_IDEAL_GAS  PV = nRT for a mixture: sum of species = P/(RT).
%   C [mol/m^3] N x 4, P_mmHg [mmHg] N x 1 (or scalar). Returns max |error| [mol/m^3].
err = max(abs(sum(C, 2) - P_mmHg(:) * p.mmHg / (p.r_gas * p.t_body)));
end
