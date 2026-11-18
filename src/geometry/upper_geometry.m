function g = upper_geometry(p)
% UPPER_GEOMETRY  Gens 0-4 (trachea to segmental bronchi).
%   g = upper_geometry(p)   -> struct (see weibel_geometry for fields)
g = weibel_geometry(p, p.gens_upper);
end
