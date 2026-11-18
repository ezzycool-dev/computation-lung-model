# Methods

This document derives every equation in the model, lists each assumption, and names the
source of every constant. Bracketed tags refer to [`REFERENCES.md`](REFERENCES.md).

## 1. Domain

| Region | Volume | Treatment |
| --- | --- | --- |
| Nose/mouth → larynx (extrathoracic) | 72 mL [Brewer] | resistance $K_1$ (calibrated); gas carried in the transport chain |
| **Generations 0–4** (trachea → segmental bronchi) | 22.5 mL | **detailed: geometry, resistance, inertance, 30 transport segments** |
| Generations 5–16 | 55.5 mL | boundary: Poiseuille resistance, gas carried in the chain |
| Alveoli | — | boundary: linear compliance, fixed alveolar gas [Guyton] |

The 150 mL anatomical dead space [Guyton] is split as 72 mL extrathoracic, with the
remaining 78 mL shared between Gens 0–4 and Gens 5–16 in proportion to their Weibel
volumes (50.5 vs 124.9 mL). Gens 0–4 therefore hold 78 × 0.288 = 22.5 mL.

## 2. Geometry

Weibel's symmetric model A [Weibel; Murayama]: generation $k$ contains $n_k = 2^k$ identical
tubes. With diameter $d_k$ and length $L_k$:

$$A_k = n_k \pi d_k^2 / 4, \qquad V_k = A_k L_k$$

| Gen | $n$ | $d$ (cm) | $L$ (cm) | $A$ (cm²) | $V$ (mL) |
| --: | --: | --: | --: | --: | --: |
| 0 | 1 | 1.80 | 12.00 | 2.54 | 30.5 |
| 1 | 2 | 1.22 | 4.76 | 2.34 | 11.1 |
| 2 | 4 | 0.83 | 1.90 | 2.16 | 4.1 |
| 3 | 8 | 0.56 | 0.76 | 1.97 | 1.5 |
| 4 | 16 | 0.45 | 1.27 | 2.54 | 3.2 |

Total cross-section barely changes across Gens 0–4, then grows about 70-fold by Gen 16. By
continuity ($u = Q/A$), air moves fastest in Gens 0–4, which is why this region carries the
largest resistance per unit volume.

## 3. Resistance

**Poiseuille baseline** (laminar, fully developed flow; [Feher] Eqn 6.2.1), with tubes in
parallel within a generation:

$$R_{P,k} = \frac{8 \mu L_k}{\pi r_k^4 n_k}$$

**Oscillatory flow (Womersley).** For flow oscillating at $\omega = 2\pi/T$ in a rigid tube,
the Navier–Stokes solution gives the impedance [Womersley]

$$Z = \frac{i \omega \rho L}{\pi r^2 \left(1 - F_{10}\right)},\qquad F_{10} = \frac{2 J_1(\Lambda)}{\Lambda J_0(\Lambda)},\qquad \Lambda = i^{3/2}\,\mathrm{Wo},\qquad \mathrm{Wo} = r\sqrt{\rho\omega/\mu}$$

$\Phi_W = \mathrm{Re}(Z)/R_P$ scales the resistance, and $L_k = \mathrm{Im}(Z)/(\omega n_k)$ is the
inertance. At a 4 s breath $\mathrm{Wo} \le 2.7$, so $\Phi_W \approx 1.0$–$1.05$.

**Entrance flow (Pedley).** Each airway is too short for the velocity profile to fully
develop after a bifurcation. Pedley, Schroter & Sudlow's correction [Pedley] is

$$\Phi_P = 0.327\sqrt{\mathrm{Re}\, d/L},\qquad \mathrm{Re} = \rho u d/\mu$$

**Combined:** $R_k(Q) = R_{P,k}\max(\Phi_W, \Phi_P, 1)$. At peak flow $\Phi_P$ = 3.9–7.2 in
Gens 0–4, the dominant correction.

**Nose to larynx.** Rohrer form $\Delta P = K_1 Q + K_2 Q|Q|$ [Feher Eqn 6.2.10], with $K_2 = 0$
because Re < 2000 at rest. $K_1$ is the only calibrated constant. Guyton states that about
1 cmH₂O of alveolar pressure moves 0.5 L in the 2 s of quiet inspiration [Guyton]. For a
half-sine inspiration the peak flow is $Q_{pk} = \pi V_T / (2T_i)$, so

$$K_1 = \frac{\Delta P_{alv}}{Q_{pk}} - \frac{\Delta P_{tree}(Q_{pk})}{Q_{pk}} = 1.59\ \text{mmHg·s/L}$$

## 4. Coupled breathing

Mass and momentum for the whole airway, with the lung as the downstream boundary:

$$\frac{dV_A}{dt} = Q$$

$$L_{tot}\frac{dQ}{dt} = P_{atm} - P_{alv} - K_1 Q - \sum_{k=0}^{4} R_k(Q)\,Q - R_{5\text{–}16}\,Q$$

$$P_{alv} = P_{atm} + P_{pl}(t) - P_{pl,FRC} + \frac{V_A - V_{FRC}}{C_L}$$

The pleural drive is a half-cosine rise over 2 s and fall over 2 s from −5 cmH₂O [Guyton].
Its amplitude is rescaled until $V_T = 500$ mL (2–3 passes). The system is stiff (airway
$L/R$ ~ ms vs a 4 s breath), so it is solved with `ode15s` for three breaths to reach a
periodic steady state.

## 5. Pressure profile

The momentum balance on generation $k$ (Newton's second law on the gas column):

$$P_k - P_{k+1} = R_k(Q)\,Q + L_k\,\frac{dQ}{dt}$$

Starting from $P_0 = P_{atm} - K_1 Q$ at the trachea entry, this gives the absolute pressure at
every branch point. Check 3 confirms the path closes on the alveolar pressure to 1e-13 mmHg.

## 6. Gas transport

The airway is rigid, so the same volumetric flow $Q$ passes every cross-section. Gas
concentrations change only by convection, with no reaction and no wall exchange (no alveoli
before Gen 17 [Feher]).

The chain of $N$ well-mixed segments runs from the mouth to the alveoli: 96 extrathoracic,
30 in Gens 0–4 and 74 in Gens 5–16, all of 0.75 mL. Each gas obeys an upwind (donor-cell)
balance:

$$V_j \frac{dC_{s,j}}{dt} = |Q|\left(C_{s,\text{up}} - C_{s,j}\right),\qquad
C_{s,\text{up}} = \begin{cases} C_{s,j-1} & Q \ge 0 \\ C_{s,j+1} & Q < 0 \end{cases}$$

Inspired gas is humidified tracheal air (P$_{O_2}$ 149.4, P$_{CO_2}$ 0.2, P$_{N_2}$ 563.3,
P$_{H_2O}$ 47 mmHg [Feher]). On expiration the deepest segment receives alveolar gas
(100 / 40 / 573 / 47 mmHg [Guyton]).

As the segment volume → 0 the chain approaches the advection equation
$\partial C/\partial t + Q\,\partial C/\partial V = 0$, a sharp front, which is plug flow. With
finite segments the front has a small numerical spread. Check 12 shows the dead space is
unchanged (< 0.01 mL) when the segment count doubles.

**Compartment form.** Summing over the 30 segments of Gens 0–4 gives the report's equation

$$\frac{dN_s}{dt} = Q_{in} C_{s,in} - Q_{out} C_{s,out}$$

## 7. Derived quantities

**Gas crossing a boundary.** The partial volume of gas $s$ in one breath is
$\int |Q|\,C_s / \sum_s C_s \, dt$, split into inspired and expired parts.

**Bohr dead space** [Robertson]: $V_D = V_E\,(P_{A,CO_2} - \bar P_{E,CO_2})/(P_{A,CO_2} - P_{I,CO_2})$,
with $\bar P_{E,CO_2}$ the flow-weighted mean expired CO₂ at the measuring point. At the mouth
it recovers the full 150 mL. At the trachea it recovers the 78 mL intrathoracic chain.

**Fick check.** With fixed alveolar gas, net uptake per breath must equal
$(V_T - V_D)(P_I - P_A)/P_B$: 22.75 mL O₂ and 18.33 mL CO₂. The model gives 22.76 and 18.33 mL.

## 8. Constants

| Symbol | Value | Source |
| --- | --- | --- |
| $P_{atm}$ | 760 mmHg | [Guyton] |
| $\rho$, $\mu$ (air, 37 °C) | 1.13 kg/m³, 19.1 µPa·s | [Feher] |
| $V_T$, breath period | 500 mL, 4 s (2 s in) | [Guyton] |
| Dead space | 150 mL total, 72 mL extrathoracic | [Guyton], [Brewer] |
| $C_L$ | 200 mL/cmH₂O | [Guyton] |
| $V_{FRC}$ | 2.3 L | [Guyton] |
| $P_{pl,FRC}$ | −5 cmH₂O | [Guyton] |
| Weibel $d_k$, $L_k$ | table above | [Weibel], [Murayama] |
| Segments in Gens 0–4 | 30 × 0.75 mL | convergence study (check 12) |

## 9. Limitations

- Symmetric, rigid airways. Real branching is irregular [Davies], and airway diameter grows
  with lung volume [Feher; Robertson].
- Laminar quiet breathing only.
- Fixed alveolar gas and a linear lung. These are boundary conditions here; the team model
  replaced them with dynamic gas-exchange and mechanics modules.
- No diffusion or cardiogenic mixing within segments [Robertson].
