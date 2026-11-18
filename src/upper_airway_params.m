function p = upper_airway_params()
% UPPER_AIRWAY_PARAMS  Every constant used by the upper-airway model.
%   p = upper_airway_params()
%
% One value per line, with its unit and its source in the comment.
% Sources (full references in docs/REFERENCES.md):
%   [Guyton]   Hall JE. Guyton and Hall Textbook of Medical Physiology, 12th ed. 2011.
%   [Feher]    Feher JJ. Quantitative Human Physiology, 2nd ed. 2017.
%   [Weibel]   Weibel human model A, tabulated in Murayama et al., J Drug Deliv 2019,
%              Table 1 (doi:10.1155/2019/3921426).
%   [Brewer]   Brewer LM, Orr JA. Anatomic dead space cannot be predicted by body
%              weight. Utah State University; 2007.
%   [Pedley]   Pedley TJ, Schroter RC, Sudlow MF. Respir Physiol 1970;9:387-405.
%
% Units used everywhere outside the SI physics helpers:
%   pressure mmHg (absolute unless the name says gauge) | flow L/s (+ = inspiration)
%   volume L | time s | gas concentration mol/m^3, columns O2 CO2 N2 H2O

%% Unit conversions
p.mmHg  = 133.322;              % Pa per mmHg
p.cmH2O = 98.0665;              % Pa per cmH2O
p.cm2mm = p.cmH2O / p.mmHg;     % mmHg per cmH2O (0.7356)

%% Air at body temperature
p.p_atm   = 760;                % mmHg   sea-level barometric pressure [Guyton]
p.rho_air = 1.13;               % kg/m^3 air density at 37 C [Feher]
p.mu_air  = 19.1e-6;            % Pa s   air viscosity at 37 C [Feher]
p.r_gas   = 8.314;              % J/(mol K)
p.t_body  = 310.15;             % K      37 C
p.species = {'O2', 'CO2', 'N2', 'H2O'};
% Humidified tracheal air (inspired gas once the nose has warmed and wetted it) [Feher]
p.pp_trach = [149.4 0.2 563.3 47];                        % mmHg
p.c_trach  = p.pp_trach * p.mmHg / (p.r_gas * p.t_body);  % mol/m^3

%% Breathing pattern (quiet adult breathing)
p.t_cycle      = 4;             % s    one breath (15 breaths/min)
p.t_in         = 2;             % s    inspiration ("the 2 seconds" of quiet inspiration) [Guyton]
p.t_out        = p.t_cycle - p.t_in;
p.n_time       = 201;           % points per breath on the output grid
p.n_breaths    = 3;             % breaths simulated to reach periodic steady state
p.tidal_volume = 0.500;         % L    fixed tidal volume [Guyton]
p.VT_tol       = 1e-4;          % L    tolerance on the solved tidal volume (0.1 mL)

%% Dead space
p.dead_space_vol = 0.150;       % L    anatomical dead space, ~150 mL [Guyton]
p.V_extra        = 0.072;       % L    nose/mouth to larynx (extrathoracic) [Brewer]
p.n_seg_upper    = 30;          % well-mixed segments in Gens 0-4 (0.75 mL each)
p.track_extrathoracic = true;   % carry the 72 mL nose-larynx gas in the transport chain

%% Airway geometry, Gens 0-23 [Weibel]
p.weibel_d_cm = [1.80 1.22 0.83 0.56 0.45 0.35 0.28 0.23 0.186 0.154 0.130 0.109 ...
                 0.095 0.082 0.074 0.066 0.060 0.054 0.050 0.047 0.045 0.043 0.041 0.041];
p.weibel_l_cm = [12.0 4.76 1.90 0.76 1.27 1.07 0.90 0.76 0.64 0.54 0.46 0.39 ...
                 0.33 0.27 0.23 0.20 0.165 0.141 0.117 0.099 0.083 0.076 0.059 0.050];
p.gens_upper  = 0:4;            % THIS MODEL: trachea (Gen 0) to segmental bronchi (Gen 4)
p.gens_lower  = 5:16;           % downstream conducting airways (boundary condition only)

%% Resistance
p.use_pedley  = true;           % entrance-flow correction [Pedley]
p.dP_alv_ref  = 1.0 * p.cm2mm;  % mmHg  |alveolar P| at peak quiet inspiration, ~1 cmH2O [Guyton]
p.VT_ref      = 0.5;            % L     volume that pressure moves ... [Guyton]
p.T_insp_ref  = 2;              % s     ... in this time [Guyton]
p.K2          = 0;              % mmHg/(L/s)^2 turbulent term; Re < 2000 at rest [Feher]

%% Downstream boundary: a simple lung so the airway has something to breathe into
p.C_L       = 200e-3 / p.cm2mm; % L/mmHg  lung compliance, 200 mL/cmH2O [Guyton]
p.V_FRC     = 2.300;            % L       functional residual capacity [Guyton]
p.P_pl_FRC  = -5 * p.cm2mm;     % mmHg    pleural pressure at FRC, -5 cmH2O [Guyton]
p.pp_alv    = [100 40 573 47];  % mmHg    alveolar gas returning on expiration [Guyton]

p.version = 'upper-airway-model 1.0';
end
