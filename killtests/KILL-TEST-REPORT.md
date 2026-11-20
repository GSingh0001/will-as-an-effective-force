# Kill Test Report

Harness: killtests/trajectories.py (simulator) and killtests/run_killtests.py (tests). Machine-readable output: killtests/results.json. Every number below is taken from that file, which this script writes in the same run.

Fixed random seeds: {'validation': 20261118, 'validation_eta05': 20261119, 'K1': 20261120, 'K2': 20261121, 'K3': 20261122}. Generator: numpy default_rng (PCG64). Integration: Euler-Maruyama, dt = 0.005 in units of the inverse level splitting (hbar = 1). No model repair is attempted in this job, per the master prompt.

## Models under test

- Model A: continuous measurement of sigma_z. Conditional states follow the diffusive stochastic master equation with measurement strength gamma and efficiency eta (see trajectories.py docstring for the equation). The stochastic element is a Wiener process representing the measurement record. This is the validated reference any repaired coupling model must match in behaviour.

- Model B: the draft coupling as written, V = -W |1><1| with W = 1 constant, added to the Hamiltonian and evolved unitarily (exactly, by diagonalisation). Deterministic and Hermitian.

- Model C (diagnostic only, K3): a deterministic state-dependent population flow dp/dt = 2 p(1-p)(2p-1). It is neither A nor B. It exists to measure, in numbers, the ensemble nonlinearity that any deterministic localisation mechanism of this type produces.

## Validation of Model A against the Lindblad equation

Method: initial state with |c1|^2 = 0.3, H = (1.3/2) sigma_z, gamma = 1, T = 4, 4000 trajectories. The trajectory ensemble mean is compared with RK4 Lindblad evolution at 9 record times. The Lindblad solution is itself checked against the analytic coherence rho_01(t) = rho_01(0) exp(-i omega t) exp(-2 gamma t).

- eta = 1.0: max |ensemble mean - Lindblad| = 0.0074; Monte Carlo scale 1/sqrt(N) = 0.0158; Lindblad vs analytic = 3.42e-11. Verdict: agreement within Monte Carlo error. Pass.

- eta = 0.5: max deviation = 0.0065 against the same Monte Carlo scale 0.0158. Pass. Efficiency changes the conditional states, not the ensemble mean, as required.

## K1 Born rule

Method: H = 0, gamma = 1.0, eta = 1.0, T = 8.0 (coherence suppression exp(-2 gamma T) = exp(-16)), 2000 trajectories per value of |c1|^2. Outcome of a trajectory: final population p = rho_11 above or below 1/2. Intervals are 95 percent Wilson score intervals. Localised fraction: final p within 0.01 of 0 or 1. Model B is evolved exactly; its population never moves, so its frequency is 0 or 1 by the same threshold rule and its interval is the Wilson interval for 0 or 2000 of 2000.

| |c1|^2 | A freq | A 95% CI | A in CI | A localised | B freq | B final p | B in CI |
|---|---|---|---|---|---|---|---|
| 0.10 | 0.104 | [0.092, 0.119] | True | 1.000 | 0.000 | 0.100 | False |
| 0.25 | 0.246 | [0.228, 0.266] | True | 1.000 | 0.000 | 0.250 | False |
| 0.50 | 0.487 | [0.465, 0.509] | True | 1.000 | 1.000 | 0.500 | False |
| 0.75 | 0.749 | [0.729, 0.767] | True | 1.000 | 1.000 | 0.750 | False |
| 0.90 | 0.890 | [0.876, 0.903] | True | 1.000 | 1.000 | 0.900 | False |

Model A verdict: PASS. Every Born value inside its interval: True. Final populations localise to the poles, so the threshold rule reads out a completed selection, not an intermediate state.

Model B verdict: FAIL. The deterministic pull changes only a phase in this basis. Final population equals initial population at every row (B final p column = |c1|^2 exactly), no trajectory localises except trivially, and the threshold frequency is 0 below |c1|^2 = 1/2 and 1 above it. At |c1|^2 = 0.5 the state sits on the threshold and the rule is arbitrary. A deterministic Hermitian potential therefore does not reproduce Born frequencies at any unequal superposition. Failure margins: at |c1|^2 = 0.25 the frequency is 0.000 against 0.25; at 0.75 it is 1.000 against 0.75. These are not statistical fluctuations. They are exact.

## K2 Energy drift

Method: H = (Omega/2) sigma_x with Omega = 1, measurement along sigma_z with gamma = 1.0, eta = 1.0, started in the ground state of H (E = -0.5), 2000 trajectories plus the Lindblad ensemble. Heating rate: slope of ensemble energy E(t) = Tr(rho H) fitted over t in [0, 0.3]. Analytic ensemble value at t = 0 is gamma * Omega = 1 in model units (D[sigma_z] damps <sigma_x> at rate 2 gamma).

- Model A trajectory ensemble slope: 0.829. Lindblad slope: 0.787. Analytic instantaneous rate at t = 0: 1.000. The fitted slopes average the rate over [0, 0.3], during which the heating rate itself decays with the polarisation, so both fitted values sit below the instantaneous rate. Trajectory and Lindblad slopes agree with each other to about 5 percent, which is the cross-check. Energy then saturates toward E = 0 (the maximally mixed value for this Hamiltonian), so the qubit heats at a finite initial rate and stops at infinite temperature. Total energy is bounded because the spectrum is bounded.

- SI scale, stated assumption: level splitting hbar Omega = 1 ueV (Omega = 1.519e+09 rad/s) with gamma = Omega in model units gives an initial heating power of 1.916e-16 J/s for the single qubit. This number scales as (splitting)^2 / hbar at fixed gamma/Omega and is reported only to fix the conversion, not as a prediction for any specific device.

- Published comparison (CSL): the free-particle CSL heating rate dE/dt = 3 lambda hbar^2 m / (4 m0^2 rC^2) (Bassi, Lochan, Satin, Singh and Ulbricht, Rev. Mod. Phys. 85, 471 (2013)), with rC = 1e-7 m and m0 the nucleon mass, gives per nucleon 4.987e-44 J/s at the GRW value lambda = 1e-16 /s and 4.987e-36 J/s at Adler's value lambda = 1e-8 /s. Experimental heating and force-noise searches bound lambda from above at fixed rC, so any repaired coupling model with a free strength parameter inherits the same form of constraint: its unconditional heating must sit below the measured bounds at its own mass scaling. The toy qubit rate above is not directly comparable, because it describes monitored dephasing of one qubit, not spontaneous mass-proportional noise. The comparison fixes the form a repaired model must be checked in. It does not bound Model A.

- Model B: total energy under H + V is conserved exactly (drift 0.00e+00 over the run, zero to integration precision, and exact by construction since the evolution is a single diagonalisation). The draft's description of the coupling as an energetic sink is therefore not a property of the equation as written: a time-independent Hermitian term conserves the total energy it helps define. If W(t) is ramped externally, the energy change is work done by the ramp, set by the ramp, not by the state. Verdict: Model B neither heats nor drains. It fails K2 as a sink and passes K2 as a non-heater, and both statements are the same fact: unitary dynamics does no selection work. Model A verdict on K2: heating measured, bounded, and accounted by the measurement channel. No divergence.

## K3 Ensemble linearity

Method: rho_a with |c1|^2 = 0.2, rho_b with |c1|^2 = 0.8, mixture rho_mix = (rho_a + rho_b)/2, H = (1.3/2) sigma_z, T = 3.0. Way 1: 4000 trajectories initialised half in rho_a and half in rho_b, averaged. Way 2: Lindblad evolution of rho_mix. Deviation: trace distance between the two final states. Any ensemble-level nonlinearity would permit superluminal signalling with entangled pairs (Gisin, Helv. Phys. Acta 62, 363 (1989)) and is disqualifying.

- Model A: trace distance 0.0063, Monte Carlo scale 1/sqrt(N) = 0.0158. Consistent with zero nonlinearity at the sampling level. PASS.

- Model B: trace distance 5.55e-17. Unitary evolution is exactly linear in rho. PASS on K3. Model B fails K1, not K3. The two tests separate: linearity is necessary, not sufficient.

- Model C diagnostic (deterministic state-dependent pull): trace distance 0.322 between evolve-then-mix and mix-then-evolve. This is the size of the ensemble nonlinearity that a deterministic localisation flow of this generic form produces at these parameters, four orders of magnitude above the Monte Carlo scale. A repaired model built that way would be disqualified by K3 even if its single-state trajectories localised.

## Verdicts

| Test | Model A (continuous measurement) | Model B (draft deterministic pull) |
|---|---|---|
| Validation vs Lindblad | PASS, max dev 0.0074 (eta=1), 0.0065 (eta=0.5) | not applicable |
| K1 Born rule | PASS, all five frequencies inside 95% CIs | FAIL, frequencies exactly 0 or 1, populations unmoved |
| K2 Energy | heating slope 0.787 (analytic 1), bounded, saturates | total-energy drift 0.0e+00: no sink, no heating |
| K3 Linearity | PASS, trace distance 0.0063 vs MC scale 0.0158 | PASS, trace distance 5.6e-17 (exact linearity) |

## Requirements a repaired model must meet to pass K1

Stated as mathematical requirements, in the order K1 tests them:

1. Trajectory-level selection. Individual conditional states must approach the measurement eigenstates (|rho_11 - 0| or |rho_11 - 1| small) with probability approaching 1 as gamma t grows. Ensemble dephasing alone does not pass K1, because the threshold readout then acts on unlocalised states.

2. Martingale population. The conditional population p(t) = rho_11(t) must be a martingale under the model's own noise, E[p(t)] = p(0), with absorbing behaviour at the poles. Absorption plus the martingale property forces the absorption probability at p = 1 to equal p(0) = |c1|^2, which is the Born rule. Any drift term in dp that does not vanish in conditional expectation breaks this and shows up in K1 as a frequency offset at unequal superpositions. The draft pull has the opposite structure: dp = 0 identically.

3. A genuine stochastic element. The noise must enter the conditional state equation (for diffusive models, a Wiener term with the innovation structure A rho + rho A - 2<A> rho, or a jump process with the corresponding structure). A deterministic term, Hermitian or not, that is a fixed function of the state either leaves populations unmoved in its own eigenbasis (Model B, K1 fail) or localises deterministically to the nearer pole, giving threshold frequencies 0 and 1 with a step at |c1|^2 = 1/2 (Model C form, K1 fail at every unequal point and K3 fail as well).

4. Linear ensemble. Averaging the conditional equation over the noise must give a Lindblad-form master equation, linear in rho with completely positive generators. This is the K3 requirement and is what separates stochastic collapse models (GRW, CSL, continuous measurement) from deterministic nonlinear modifications (Gisin 1989).

5. Basis from the coupling. The measured observable whose eigenstates are the absorbing poles must be the operator in the stochastic term, stated explicitly. In the draft, |1> is singled out by hand in V, which assumes the outcome basis and, with a fixed sign of W, the favoured outcome. A repaired model must get the basis from a stated system-environment or system-apparatus interaction, not from the answer.

6. Bounded, accounted energy. The unconditional energy drift implied by the Lindblad ensemble must be computed and placed against measured heating bounds at the model's mass scaling (K2). For diffusive measurement of an observable that does not commute with H, a nonzero heating rate is generic (measured above: slope 0.787 at gamma = Omega = 1). A model that claims an energy sink must show the sink in the ensemble equation. A Hermitian potential cannot.

## Open items in this job

- The simulator covers one qubit, Markovian diffusive measurement, and efficiency-independent ensemble checks at eta = 1 and 0.5. Jump (photon-counting) unravellings, non-Markovian records, and multipartite systems are not implemented. K1 to K3 as specified do not require them.

- K2's SI conversion assumes a 1 ueV splitting to fix units. The CSL figures quoted use rC = 1e-7 m and the GRW and Adler lambda values. Both are stated inputs. Current experimental upper bounds on lambda at this rC are stronger than Adler's value in several channels. This report does not quote a single bound number because the bound is a curve in the lambda-rC plane, and reducing it to one number without the curve would misstate it. Left OPEN as a deliberate omission, with the source (Bassi et al. 2013 and subsequent bound compilations) named.
