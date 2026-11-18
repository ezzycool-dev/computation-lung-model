function dP = upper_dP(Q_V, p, g)
% UPPER_DP  Resistive pressure drop across Gens 0-4.
%   Q_V [L/s] scalar or vector (signed, + = inspiration)
%   dP  [mmHg], same size and sign as Q_V
% Law: Hagen-Poiseuille resistance scaled by max(Womersley, Pedley, 1);
% generations in series add, tubes within a generation are in parallel.
dP = reshape(sum(generation_drops(Q_V, p, g), 2), size(Q_V));
end
