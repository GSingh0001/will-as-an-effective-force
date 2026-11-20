# Born, Energy and No-Signalling Tests of Observational Coupling Models: A Reproducible Falsification Study

Ganga Singh
Independent Researcher
ORCID: 0009-0004-2191-356X

## Abstract

Candidate models in which an observational or coupling term is added to quantum dynamics are usually assessed first by interpretation: what the new term might mean, and what it might explain. This paper argues for the opposite order. Before interpretation, a candidate should be subjected to three cheap, decisive, fully numerical tests that any such model must pass to remain a candidate at all. The tests are: (K1) recovery of Born rule outcome frequencies from individual conditional trajectories; (K2) an energy audit, requiring that the model's unconditional energy behaviour be computed, bounded, and consistent with the role claimed for it; and (K3) a finite ensemble-consistency test, comparing two orders of mixing and evolving for one specified decomposition, motivated by the ensemble-linearity condition associated with Gisin's no-signalling argument. We implement the three tests as a small frozen harness for a single qubit, using Euler-Maruyama trajectory integration at fixed random seeds, and we publish the harness protocol in enough detail that an outside group can rerun it from this paper alone. We then apply the harness to a reference model and to a family of candidate couplings. Standard continuous measurement of a qubit observable (Model A) passes all three tests. A deterministic Hermitian coupling of the form V = -W |1><1| (Model B) fails K1 exactly: its outcome frequencies are 0, 0, 1, 1, 1 against Born probabilities 0.10, 0.25, 0.50, 0.75, 0.90, because the coupling changes only a phase and no population ever moves. The same model conserves energy exactly, so it also fails in the role of an energy sink that its motivating description assigns to it. A deterministic state-dependent pull, included as a diagnostic (Model C), fails the tested K3 ensemble-consistency condition with a trace distance of 0.3223 between the two orders of mixing and evolving, against a Monte Carlo scale of 0.0158. Four further frozen-run candidates (a drift-biased diffusion, a naive hit model, and two coloured-noise constructions) each fail K1 at most or all of the five test points, while a Tilloy-Diosi form passes the form checks applied to it. Finally, we show that a stochastic repair of Model B, obtained by replacing the deterministic term with a diffusive innovation structure satisfying the requirements stated in Section 5.1, is not a new model at all: after the parameter substitution stated in the Appendix, its stochastic master equation is identical to standard continuous measurement theory in the form given by Wiseman and Milburn, equivalently the Belavkin filtering equation. The harness therefore does double duty. It falsifies candidate couplings cheaply, and when a candidate is repaired until it passes, it can reveal that the repair has quietly reproduced known theory. A numerical pass in K1 or K3 is a compatibility statement at the finite sensitivity of the run, not a proof of exact invariance, except where a verdict follows directly from the model's equations. All runs reported here were re-executed from frozen code in a separate verification pass, and every verdict reproduced exactly.

## 1 Introduction

New terms get added to quantum dynamics for many reasons: as fundamental modifications of the Schroedinger equation, in the tradition of spontaneous collapse models; as effective descriptions of measurement, feedback, or monitoring; or as couplings between a quantum system and an external agent, record, or observer-like degree of freedom. Whatever the motivation, the resulting equations share a common burden: they must reproduce the statistical content of ordinary quantum mechanics where it has been tested, they must not inject or remove energy in ways that contradict their own stated purpose or experimental heating bounds, and they must not permit superluminal signalling. These are not interpretive demands. They are arithmetic demands, and they can be checked arithmetically.

In practice, candidate models are often discussed at length before anyone checks the arithmetic: the term is written down, its meaning debated, its consequences explored, and only later does someone simulate it and discover that its trajectories do not localise, or that its outcome frequencies are a step function rather than the Born rule, or that its ensemble dynamics would allow signalling. This is an expensive order of operations, because interpretive debate is unbounded while the arithmetic checks take minutes. They should come first.

This paper makes four contributions. First, it defines the three tests precisely enough to be implemented independently. Second, it describes a frozen reference implementation for a single qubit, with fixed seeds and acceptance criteria, and a rerun protocol an outside group can follow from the paper alone. Third, it reports verdicts for a reference model, for a deterministic coupling that fails exactly and instructively, for a diagnostic flow that quantifies how badly deterministic localisation violates linearity, and for five further candidates graded in frozen runs. Fourth, it reports an identity result: a repair that makes the failed deterministic coupling pass the tests has a stochastic master equation that is term-by-term identical to standard continuous measurement theory after an explicit parameter substitution. The identity converts a negative result (this coupling fails) into a closure result for that repaired equation: the diffusive repair considered here contains nothing new. A researcher who arrives at the same equations by a different route has rediscovered continuous measurement, and is better off knowing that at the start of the project than at the end.

Two scope remarks. First, the harness is deliberately small: one qubit, Markovian diffusive dynamics, and two detection efficiencies suffice for every failure reported here, because the failures are structural rather than parametric, and several are exact. Second, this paper makes no claim about which modification of quantum dynamics, if any, is correct. It is a paper about a sieve. Models that fail it are finished, whatever their motivation; models that pass it have earned further discussion, and nothing more.

## 2 The three tests

Each test is stated with its protocol, its acceptance criterion, and the reason the criterion is the right one. The numerical settings quoted here are those of the frozen harness described in Section 3; the criteria themselves are stated so that other implementations can adopt them at other settings.

### 2.1 K1: Born rule statistics from trajectories

Protocol. Set the Hamiltonian to zero in the measured basis, so that no unitary dynamics competes with the coupling. Initialise pure states with excited-state population p0 = |c1|^2 at the five values 0.10, 0.25, 0.50, 0.75, 0.90. Evolve a large ensemble of conditional trajectories (2000 per value in the frozen harness) long enough that the measurement-induced coherence suppression is effectively complete (in the harness, gamma = 1 and T = 8, so coherences are suppressed by exp(-16)). For each trajectory, record the final population p and classify the outcome as 1 if p > 1/2 and 0 otherwise. Also record the localised fraction, the proportion of trajectories whose final p lies within 0.01 of 0 or 1. Compare the empirical frequency at each p0 with the Born probability p0 using a 95 percent Wilson score interval.

Acceptance. Within the sensitivity of 2000 trajectories per point, a model passes K1 if the Born value lies inside the 95 percent interval at all five test points, and if the localised fraction is effectively 1, so that the threshold rule is reading completed selections rather than intermediate states.

Scope of a numerical pass. The finite trajectory ensemble does not establish exact Born invariance. It detects departures large enough to move the observed frequencies outside the stated intervals at the chosen trajectory count, and it can miss smaller departures. For Model A and for the repair in Section 5, Born absorption also follows from the martingale structure built into the conditional population equation. The numerical K1 run is therefore a code-consistency test of that construction and of its implementation, not an independent empirical proof of the Born rule. Model B is different: its K1 failure is exact because its equations leave every population unchanged.

Why trajectories, and why this is decisive. Ensemble dephasing is not enough: a model whose average state loses coherence shows a diagonal ensemble density matrix, but if its individual conditional states do not approach the eigenstates, no single run ever produces a definite outcome, and the threshold readout is an artefact of the analyst rather than a prediction of the model. The mathematical content of the test is the martingale property. If the conditional population p(t) is a martingale under the model's own noise, with absorption at the poles p = 0 and p = 1, then the absorption probability at p = 1 equals p(0), which is the Born rule. Any drift in p that does not vanish in conditional expectation breaks the frequency match at unequal superpositions by an amount that grows with trajectory count rather than shrinking. A K1 failure far outside the intervals is therefore evidence of a structural defect, while a pass leaves open departures below the finite sampling sensitivity.

### 2.2 K2: energy behaviour

Protocol. Choose a Hamiltonian that does not commute with the measured observable, so that measurement backaction does work on the system. In the harness, H = (Omega/2) sigma_x with Omega = 1, measurement is along sigma_z with gamma = 1 and efficiency eta = 1, and the system starts in the ground state of H, with energy E = -0.5. Evolve 2000 trajectories and the corresponding ensemble master equation. Fit the slope of the ensemble energy E(t) = Tr(rho H) over an early window (t in [0, 0.3] in the harness) to obtain a heating rate, and track the energy over the full run to establish whether it is bounded. For a deterministic candidate, compute the total energy drift under the full Hamiltonian including the coupling term.

Acceptance. K2 is an audit rather than a single threshold. A model passes if (i) its unconditional energy behaviour is computed and reported, (ii) that behaviour is bounded, or its growth is explicitly characterised, and (iii) the behaviour matches the role claimed for the coupling. A model advertised as an energy sink must show draining in its ensemble equation. A model that heats must show a finite rate that can be placed against experimental heating bounds at the model's mass scaling. Exact conservation is a legitimate outcome for a unitary model, but it disqualifies any claim that the coupling performs selection work energetically, because a time-independent Hermitian term conserves the total energy it helps define. If an external agent ramps the coupling in time, the energy change is work done by the ramp, set by the ramp and not by the state, and it must be accounted as such.

Scope of the audit. An energy increase has the force of a rejection only against the model's stated accounting assumptions. It rejects a model that promises conservation for an isolated system, or that claims a specific energetic role contradicted by its ensemble equation. It does not by itself reject an effective open-system model, for which exchange with an apparatus or environment may be part of the model. The converse also holds: bounded or otherwise acceptable mean energy behaviour in this test does not prove a complete conservation ledger. K2 reports the unconditional mean energy for the stated Hamiltonian and settings; it does not, by itself, account for every work, heat, record, or entropy flow in an implementation.

Why an audit. Collapse and measurement models generically heat: for diffusive measurement of an observable that does not commute with H, a nonzero heating rate is the rule, and in the harness the analytic instantaneous rate at t = 0 is gamma times Omega, equal to 1 in model units. Spontaneous collapse models inherit the same structure at their own scaling: the free-particle heating rate of continuous spontaneous localisation (CSL) is dE/dt = 3 lambda hbar^2 m / (4 m0^2 rC^2), which at rC = 1e-7 m gives 4.987e-44 J/s per nucleon at the GRW value lambda = 1e-16 /s and 4.987e-36 J/s at Adler's value lambda = 1e-8 /s, following Bassi et al. (2013). (The identity run recomputed the same formula and obtained 4.99e-44 and 4.99e-36 J/s, the same values at its printed precision.) Experimental searches bound such rates from above, so any model with a free strength parameter must eventually present its heating in this form. K2 forces that presentation at the outset, at toy scale, where it costs nothing.

### 2.3 K3: tested ensemble consistency (the finite Gisin condition)

Protocol. Take two initial states, rho_a with |c1|^2 = 0.2 and rho_b with |c1|^2 = 0.8, and their equal mixture rho_mix = (rho_a + rho_b)/2. Evolve the mixture in two ways. Way 1: initialise trajectories half in rho_a and half in rho_b (4000 trajectories in the harness), evolve each under the model's conditional dynamics, and average. Way 2: evolve rho_mix directly under the model's ensemble equation. Compute the trace distance between the two final states, and compare it with the Monte Carlo scale 1/sqrt(N), which is 0.0158 at N = 4000.

Acceptance and scope. A model passes this K3 test if the trace distance for this specified decomposition is consistent with the Monte Carlo scale, that is, if any ensemble nonlinearity for that decomposition is below the sampling floor of the test. A trace distance far above the scale is a measured nonlinearity on the tested decomposition, not a numerical error. A pass is evidence of ensemble consistency on the decomposition tested, at that sampling floor. It is not a proof of linearity over all density operators, all ensemble decompositions, or entangled extensions. For Model B, a stronger statement is available analytically because unitary conjugation is linear in rho; the tiny numerical value checks that implementation rather than supplying the reason for linearity.

Why linearity is disqualifying when it fails. Gisin showed that deterministic nonlinear modifications of the Schroedinger equation permit superluminal signalling with entangled pairs: the reduced dynamics of one member comes to depend on how the distant member's state was prepared, because different ensemble decompositions of the same density operator evolve differently (Gisin 1989). The two-ways construction measures decomposition dependence for the particular pair of states and mixture used here. That is the finite condition tested by K3. It should not be paraphrased as a proof of the full no-signalling condition for arbitrary states, decompositions, or composite systems. A large failure is nevertheless serious: for the deterministic constructions considered here, it exhibits the decomposition dependence that motivates the Gisin objection. Note the asymmetry the results below make vivid: linearity is necessary, not sufficient. A model can pass K3 exactly and still fail K1 completely.

## 3 The harness

### 3.1 Implementation

The reference implementation consists of a trajectory simulator and a test runner. Conditional states for the diffusive models follow the stochastic master equation for continuous measurement of sigma_z with strength gamma and efficiency eta, integrated by the Euler-Maruyama method with time step dt = 0.005 in units of the inverse level splitting, with hbar = 1. Deterministic candidates given by a time-independent Hamiltonian term are evolved exactly, by diagonalisation, so that no integration error can be blamed for their behaviour. The ensemble master equation used in K2 and K3 is integrated independently (RK4 in the validation runs) and is also checked against closed-form solutions where they exist.

Randomness is fixed and reported. The generator is the numpy default generator (PCG64). The frozen seeds are: validation 20261118, validation at eta = 0.5 20261119, K1 20261120, K2 20261121, K3 20261122. Trajectory counts are 4000 for validation and K3, and 2000 per test point for K1 and K2. Acceptance criteria, thresholds, and trajectory counts were fixed before the candidate models were graded, and no model repair was attempted inside the harness job: a candidate is graded as written, and repair is a separate, later question (taken up in Section 5).

### 3.2 Validation of the reference model

Before any candidate is graded, the harness validates its reference dynamics. Model A (Section 4.1) is simulated from an initial state with |c1|^2 = 0.3, with H = (1.3/2) sigma_z, gamma = 1, T = 4, and 4000 trajectories. The ensemble mean is compared with RK4 integration of the Lindblad equation at nine record times, and the Lindblad solution is itself checked against the analytic coherence rho_01(t) = rho_01(0) exp(-i omega t) exp(-2 gamma t). At eta = 1.0 the maximum deviation from the Lindblad solution is 0.0074, against a Monte Carlo scale of 0.0158; the Lindblad solution matches the analytic coherence to 3.42e-11. At eta = 0.5 the maximum deviation is 0.0065. Efficiency changes the conditional states and not the ensemble mean, as it must. The harness thus measures what it claims to measure, at a known sampling scale, before it is allowed to fail anyone.

### 3.3 Rerunning the harness from this paper alone

An outside group can reproduce every verdict in this paper as follows. (1) Implement the diffusive stochastic master equation for a monitored sigma_z channel with parameters gamma and eta, or configure any standard quantum trajectories code equivalently. (2) Fix the protocols of Section 2 exactly: the five K1 populations, the K2 Hamiltonian and fitting window, the K3 states and mixture weights. (3) Use trajectory counts at or above those of Section 3.1 and any seed; the verdicts are structural and seed-independent, as Section 4.6 demonstrates: an independent re-execution at different seeds (20270110 to 20270119) on a different machine reproduced every verdict, with frequency lists differing by at most 0.0155 at 2000 trajectories per point. (4) Implement deterministic candidates exactly, by diagonalisation. (5) Apply the acceptance criteria of Section 2 as written. The regrade of Section 4.6 was performed from byte-identical copies of the frozen code, and a second execution of the repair runner reproduced its results file with no differing byte.

One environment note from the regrade: the re-execution ran under numpy 1.26.4 while the original runs used numpy 2.0 or later, and one candidate called a function renamed between those versions (np.trapezoid in numpy 2, np.trapz before it). The regrade supplied a runtime alias and did not modify the candidate code. All verdict-bearing frequencies were exactly identical across the runs; remaining differences were confined to analytic floating-point quantities at or below 3.1e-16, plus one rig-validation quantity differing by 2.8e-17.

## 4 Models and verdicts

### 4.1 Model A: standard continuous measurement (reference)

Model A is the diffusive stochastic master equation for continuous measurement of sigma_z, with measurement strength gamma and efficiency eta. It is the reference any repaired coupling model must match in behaviour. Its K1 results at gamma = 1, eta = 1, T = 8, 2000 trajectories per point, are frequencies 0.104, 0.246, 0.487, 0.749, 0.890 at Born values 0.10, 0.25, 0.50, 0.75, 0.90, with 95 percent Wilson intervals [0.092, 0.119], [0.228, 0.266], [0.465, 0.509], [0.729, 0.767], [0.876, 0.903]. Every Born value lies inside its interval, and the localised fraction is 1.000 at every point: trajectories complete their selection, and the threshold rule reads the completed result. This is a finite-sample consistency result at 2000 trajectories per point, not a demonstration of exact Born invariance. K1: PASS.

On K2, Model A's trajectory ensemble gives a fitted heating slope of 0.829 and the Lindblad ensemble 0.787, against an analytic instantaneous rate of 1.000 at t = 0. The fitted slopes average the rate over the fitting window, during which the rate decays with the polarisation, so both sit below the instantaneous rate; they agree with each other to about 5 percent, which is the cross-check the test is designed to produce. The energy then saturates toward E = 0, the maximally mixed value: the qubit heats at a finite initial rate and stops at infinite temperature, and the total change is bounded because the spectrum is bounded. For scale only, the harness report converts the toy rates to SI units under a stated assumption: a 1 micro-eV splitting with gamma equal to Omega. On that assumption the analytic instantaneous rate of 1.000 converts to an initial heating power of 2.43e-16 J/s for the single qubit, while the fitted ensemble slope of 0.787 converts to 1.916e-16 J/s; the figure of 1.916e-16 J/s printed in the harness report is the conversion of the fitted slope, not of the initial rate. Both figures scale as (splitting)^2 / hbar at fixed gamma/Omega. K2: PASS as an audit.

On K3, Model A's trace distance between the two orders of mixing and evolving is 0.0063, against a Monte Carlo scale of 0.0158: consistent with zero nonlinearity for the tested decomposition at the sampling level. No conclusion is drawn for untested decompositions or entangled extensions. K3: PASS on the tested decomposition.

### 4.2 Model B: the deterministic coupling V = -W |1><1|

Model B is a candidate observational coupling in its draft form: a term V = -W |1><1| with W = 1 constant, added to the Hamiltonian and evolved unitarily. It is deterministic and Hermitian. In the harness it is evolved exactly.

Model B fails K1 exactly. In the sigma_z basis the term is diagonal, so it changes only a phase; the population of every state is constant in time. The final population equals the initial population at every test point, no trajectory localises except trivially, and the threshold frequency is 0 for every p0 below 1/2 and 1 for every p0 above it: frequencies 0, 0, 1, 1, 1 at Born values 0.10, 0.25, 0.50, 0.75, 0.90, as the re-executed control run confirms. At p0 = 0.5 the state sits exactly on the threshold and the readout rule is arbitrary. The failure margins are not statistical: at p0 = 0.25 the frequency is 0.000 against a Born value of 0.25; at p0 = 0.75 it is 1.000 against 0.75. A deterministic Hermitian potential does not reproduce Born frequencies at any unequal superposition, and in this basis it does not even move the quantity the readout measures. K1: FAIL, exactly.

Model B's K2 result is the same fact viewed from the other side. The total energy under H + V is conserved exactly: the measured drift over the run is 0.00e+00, zero to integration precision and exact by construction, since the evolution is a single diagonalisation. The motivating description of this coupling cast it as an energetic sink, a term that would draw energy from the system as selection occurs. That description is not a property of the equation as written. A time-independent Hermitian term conserves the total energy it helps define; it neither heats nor drains. Model B therefore fails K2 in the role of a sink, and passes K2 only in the empty sense of being a non-heater. Both statements are one statement: unitary dynamics does no selection work.

On K3, Model B passes trivially and exactly: unitary evolution is linear in rho, and the measured trace distance is 5.55e-17 (5.6e-17 at the precision printed in the verdict table of the harness report). K3: PASS, exactly for this unitary model. The juxtaposition is the lesson of this model. Linearity is necessary and nowhere near sufficient; a model can be perfectly causal and perfectly incapable of producing an outcome.

### 4.3 Model C: deterministic pull (diagnostic)

Model C is neither the reference nor the candidate. It is a diagnostic constructed to put a number on the ensemble nonlinearity that deterministic localisation produces. Its dynamics is a deterministic state-dependent population flow, dp/dt = 2 p (1-p) (2p-1), which pulls every population toward the nearer pole. Such a flow does localise single states, which is what makes it a useful diagnostic: if a repair of Model B were built in this style, its trajectories would look right while its ensemble quietly became decomposition-dependent.

Run through K3 at the standard settings, Model C gives a trace distance of 0.3223 between evolve-then-mix and mix-then-evolve. (The original harness report prints this quantity as 0.322; the identity run and the independent regrade both give 0.3223, and we report the more precise value while noting the earlier rounding.) Against a Monte Carlo scale of 0.0158, the measured nonlinearity is larger by a factor of about twenty in the trace distance itself, and by a factor of about four hundred in the squared quantities. (The harness report's own wording describes the trace distance as four orders of magnitude above the sampling scale; that overstates the ratio, and we record the error rather than repeat it.) Either way, the conclusion is not marginal: this deterministic pull fails the tested K3 ensemble-consistency condition even if its single-state trajectories localise beautifully. K3: FAIL.

### 4.4 Frozen-run candidates

Four further candidates were graded in a separate frozen run series, under the same tests, and every verdict was reproduced exactly in the independent regrade of Section 4.6. In each case the five numbers are the simulated outcome frequencies at Born values 0.10, 0.25, 0.50, 0.75, 0.90.

A drift-biased diffusion, in which the conditional population dynamics carries a systematic drift in addition to its noise, gives frequencies 0.1460, 0.3235, 0.5740, 0.8180, 0.9305. It fails K1 at all five points. This is the martingale argument of Section 2.1 made numerical: a drift that does not vanish in conditional expectation displaces every frequency in the direction of the bias, by amounts far outside sampling error.

A naive hit model, in which discrete localisation events strike the state without the innovation structure of a measurement record, gives frequencies 0.2760, 0.3510, 0.5425, 0.6400, 0.7305. It fails K1 at all five points. Its K2 fitted slope is 0.1058, and its K3 trace distance is 0.0437, which is 2.8 times the Monte Carlo scale, so it shows a measurable ensemble nonlinearity on the tested decomposition as well.

A coloured-noise model with Ornstein-Uhlenbeck noise gives frequencies 0.4440, 0.4650, 0.4875, 0.5510, 0.5825. It fails K1: the frequencies are compressed toward 1/2, the signature of a coupling whose effective action averages away before selection completes. A supplementary point at correlation time tau_c = 0.1 gives 0.4740, reproducing the same compression.

A coloured-noise model in the Wong-Zakai form, where the smooth noise limit generates an additional drift, gives frequencies 0.0300, 0.1215, 0.5030, 0.8700, 0.9670. It fails K1 at four of the five points; the midpoint, where the Born value 0.50 sits near the model's own symmetry point, is the exception. Its K3 trace distance is 0.0055, a pass on the tested decomposition: like Model B, it is consistent with ensemble linearity on that test and wrong at trajectory level, from the opposite direction, overshooting the Born frequencies instead of flattening them.

A fifth candidate, a Tilloy-Diosi form tying the coupling to a gravitational self-energy structure, was assessed by form checks rather than by trajectory grading. A form check, as the term is used here and in the abstract, is a structural comparison of the candidate's equations against the Diosi-Penrose reference run, covering the ladder ratios and the form of the ensemble equation; it is not a grading of trajectories against K1 to K3. The candidate passes the form checks applied to it: its K1 form check passes, its K3 form check gives a trace distance of 0.0063, and its Diósi-Penrose ladder ratios against the reference run equal 1.0000 at every object tested, with a dust-grain Diósi-Penrose rate of 1.020e+05 /s and a ratio of that rate to the CMB reference rate of 2.846e+06. Section 5.2 returns to it: its passing form is no coincidence, because the form is identical to the Diósi-Penrose structure it is named for.

### 4.5 Verdict table

| Model | K1 Born rule | K2 Energy | K3 tested ensemble consistency |
|---|---|---|---|
| A: continuous measurement | PASS (all five Born values in 95% CIs) | Heating slope 0.787 (ensemble), bounded, saturates | PASS on tested decomposition (0.0063 vs scale 0.0158) |
| B: deterministic V = -W |1><1| | FAIL, exact (frequencies 0, 0, 1, 1, 1) | Drift 0.00e+00: no sink, no heating; fails as a sink | PASS (5.55e-17, exact linearity) |
| C: deterministic pull (diagnostic) | Not graded (diagnostic) | Not graded | FAIL (0.3223 vs scale 0.0158) |
| Drift-biased diffusion | FAIL at all five points | (see frozen record) | (see frozen record) |
| Naive hits | FAIL at all five points | Slope 0.1058 | 0.0437 (2.8 times scale) |
| Coloured noise, Ornstein-Uhlenbeck | FAIL (frequencies compressed to 1/2) | (see frozen record) | (see frozen record) |
| Coloured noise, Wong-Zakai | FAIL at four of five points | (see frozen record) | PASS on tested decomposition (0.0055) |
| Tilloy-Diosi form | Form check PASS | (see frozen record) | Form check PASS (0.0063) |

The pattern across the table is the empirical content of this paper. Every deterministic construction fails, each in a characteristic way: by not moving populations at all, by localising at the price of ensemble nonlinearity, or by drifting frequencies off the Born values. Every coloured-noise construction fails K1 as well, in one direction or the other. The constructions that pass are exactly those already known to pass, because they are the standard theory or a named relative of it. Section 5 makes that precise for the most instructive case.

### 4.6 Independent regrade

All of the above was re-executed from frozen code in a separate verification pass, on a different machine, before this paper was written. The harness was rerun as a control: Model A's K1 frequencies reproduced as 0.104, 0.246, 0.487, 0.749, 0.890, Model B's as 0, 0, 1, 1, 1, and Model C's K3 trace distance as 0.3223 against the scale 0.0158. The four candidates of Section 4.4 reproduced their frequencies exactly, digit for digit; remaining differences were confined to analytic floating-point leaves at or below 3.1e-16. The gravity-side numbers of Section 6 were likewise regenerated, with a maximum relative difference of 5.5e-16; the ladder lifetimes in that record (electron 7.3529e+23 s, C60 3.8688e+14 s, levitated sphere 2.1117e-02 s, milligram grain 9.8018e-16 s, cat 9.7247e-27 s) were identical. No verdict changed under re-execution: the verdicts are properties of the models and the tests, not of a run, machine, or seed.

## 5 The identity result

### 5.1 A stochastic repair of Model B is continuous measurement

Section 4.2 leaves a natural question. Model B fails because it has no stochastic element: its population dynamics is dp = 0 identically, the opposite of the martingale-with-absorption structure K1 requires. What happens if the deterministic term is replaced by a stochastic structure satisfying the requirements the tests impose? The harness report states those requirements in the order the tests check them: trajectory-level selection; a martingale population with absorbing poles; a genuine stochastic element with the innovation structure A rho + rho A - 2<A> rho; a linear Lindblad ensemble; a measured basis that comes from the coupling itself rather than being inserted by hand; and bounded, accounted energy. These are requirements for the repair constructed here. They are not a definition of minimality, and no admissible class or optimisation criterion is specified in the source record. Accordingly, this paper claims only that a repair meeting them was constructed, not that it is the least or unique repair.

Write down the ansatz that meets these requirements for one qubit monitored in the sigma_z basis. With hbar = 1, measurement strength gamma, and efficiency eta, it is:

d rho_c = -i [H, rho_c] dt + gamma D[sigma_z] rho_c dt + sqrt(eta gamma) H[sigma_z] rho_c dW_t,

where D[A] rho = A rho A^dagger - (1/2) {A^dagger A, rho} and H[A] rho = A rho + rho A^dagger - Tr[(A + A^dagger) rho] rho. Averaging over the noise gives the ensemble equation d rho/dt = -i [H, rho] + gamma D[sigma_z] rho. At eta = 1 the pure-state form is d|psi_c> = [-iH - (gamma/2)(sigma_z - <sigma_z>)^2] |psi_c> dt + sqrt(gamma) (sigma_z - <sigma_z>) |psi_c> dW_t. When [H, sigma_z] = 0, the population obeys dp = -4 sqrt(eta gamma) p (1-p) dW_t, a driftless diffusion with absorbing poles: the martingale K1 demands. The ensemble heating functional is dE/dt = gamma Tr[(sigma_z H sigma_z - H) rho], the quantity K2 audits.

This ansatz was graded against the harness in a dedicated run, at new seeds (20270110 to 20270119). It passes. Its validation deviation from the Lindblad solution is 0.0060 at eta = 1 and 0.0076 at eta = 0.5, against the Monte Carlo scale 0.0158. Its K1 frequencies are 0.1075, 0.2435, 0.5025, 0.7520, 0.9005, with localised fraction 1.000 at every point; the maximum absolute difference from Model A's original frequency list is 0.0155 at 2000 trajectories per point, which is sampling variation between two runs of what will turn out to be the same model. Its K2 fitted slopes are 0.779 (trajectory ensemble) and 0.787 (Lindblad ensemble), with the same analytic instantaneous rate of 1.000. Its K3 trace distance is 0.0007 against the scale 0.0158, and the deterministic diagnostic run alongside it gives 0.3223, matching Model C. Its heating matches the analytic rate structure exactly, with the exact initial rate 1.000 stated in its derivation.

The verdict of that run, however, is not that a new model has been found. It is that the ansatz is identical to standard continuous measurement theory. The reference form, from Wiseman and Milburn's Quantum Measurement and Control, is the diffusive stochastic master equation for a monitored channel with Lindblad operator c and efficiency eta:

d rho_c = -i [H, rho_c] dt + D[c] rho_c dt + sqrt(eta) H[c] rho_c dW_t.

The substitution c := sqrt(gamma) sigma_z converts the reference form into the ansatz symbol for symbol: D[sqrt(gamma) A] = gamma D[A] for real gamma >= 0, and H[c] = sqrt(gamma) H[sigma_z] because the H superoperator is homogeneous of degree one in its argument at Hermitian argument. Every term maps, and the mapping was verified term by term: the Hamiltonian drift, the dissipator, the innovation term, the split of the channel into observed and unobserved fractions through eta, the homodyne-type diffusive record, the ensemble Lindblad equation obtained by averaging, the pure-state equation at eta = 1, the population martingale (including the recorded coefficient convention gamma_J = 4 gamma under a rescaling used in the ansatz's own derivation notes), and the heating functional. The residual freedom between the two presentations is confined to choices that change no term: the choice of H, the numerical values of gamma and eta, a record gain normalisation, the strength convention gamma against gamma_J, a sign convention, and the Ito against Stratonovich presentation after the corresponding drift conversion. In quantum probability terminology, the same equation is the Belavkin filtering equation for continuous nondemolition observation, and the diffusive single-channel case reduces to the same form.

The identity rests on the term substitution, not on the numerical agreement. Numerical agreement alone would not establish identity, since two different models can agree at five frequencies within sampling error. Term identity after that substitution is the claim, and the Appendix states its assumptions, convention, parameter mapping, and precise level. We state the consequence plainly: for this coupling form, on one qubit, under diffusive monitoring of one observable, the novelty question is closed. For the diffusive repair constructed here, the route from the failed deterministic coupling arrives at continuous measurement theory. The repaired model inherits the interpretation of the theory it is identical to; attaching a new interpretation to one of its parameters does not alter a single term, and interpretations enter no term.

### 5.2 The structural comparison with collapse models

It is worth recording why the same identity does not extend to spontaneous collapse models. In mass-proportional CSL, the ensemble equation is a double commutator with the smeared mass density operator M(x), integrated over space, with two proposed universal constants: a collapse rate lambda and a correlation length r_C. The conditional dynamics is driven by a classical noise field W_t(x), spatially correlated over r_C, and no observation record is required: the noise acts whether or not anything is measured, and there is no efficiency parameter. Two structural rows are shared with the ansatz: a linear ensemble with localising conditional dynamics, and Born frequencies for the distinguished operator. Every other row differs: the state space, the operator in the noise (a single sigma_z against a continuum family M(x)), the noise (one record Wiener process against a universal field), the parameters (an arrangement-dependent gamma against universal lambda and r_C), and the scaling (no constituent-number dependence against amplification through mass itself). Reducing either construction to the other would require maps between all of these objects, and no such map exists in the ansatz. The Tilloy-Diosi candidate of Section 4.4 sits on the other side of this divide: it passes its form checks because its structure is that of the Diósi-Penrose family, a gravitational relative of the collapse programme, and its ladder ratios of 1.0000 against the reference run say exactly that. The harness sorts candidates into three bins, and the bins are the result: constructions that fail, constructions that pass because they are the standard monitored theory, and constructions that pass because they are a named collapse model. We found no fourth bin.

A corollary concerns experimental bounds. Because the ansatz is identical to continuous measurement under the assumptions stated in the Appendix, there is no differing term for a bound to attach to: bounds on CSL's lambda constrain CSL, and importing one as a constraint on the ansatz's gamma, without an equation relating gamma to lambda and r_C, would attach a number to an undefined quantity. That step is not taken here, and should not be taken anywhere. Gamma is set by the monitoring arrangement and calibrated by observed dephasing; it is not a constant of nature awaiting a bound.

## 6 Experimental context

The harness is a toy, and this section places its results against the experimental landscape that any non-toy descendant of these models would face, using only the numbers already in the research record.

For collapse models, the binding constraints are experimental and they are strong. The record of this study cites the result of Donadi et al. (2021), an underground test of gravity-related wave function collapse obtained from spontaneous radiation emission searches of the type that have made germanium detector data one of the sharpest tools in this field. That search excludes the natural parameter-free Diosi-Penrose model and bounds the collapse rate for the mass-density-coupled, gravity-related class at lambda below 5.2e-13 /s; it is a bound on that class, and not a bound on the lambda of CSL as such. A coupling model with a free strength parameter, mass-proportional or otherwise, must present its heating and radiation signatures against bounds of this kind, as a curve in its parameter plane, and not as a single number detached from the correlation length or equivalent scale at which the bound was obtained. The K2 audit of Section 2.2 is the toy-scale rehearsal of exactly that obligation: compute the unconditional energy behaviour, state the scaling, and place it against the measured bounds.

The gravitational side of the record gives the corresponding scale for Diósi-Penrose type structures. The frozen gravity numbers include the XENONnT bound on the Diósi-Penrose parameter-free length R0, quoted as a ratio of 90.0 between the bound and the parameter-free value. (Footnote: the record preserves a discrepancy in this bound's statement, and we preserve it rather than smoothing it over. The frozen loop verdict cites the XENONnT bound as 4.5e-10 m at 95 percent confidence level; the gravity-numbers job cites 4.9e-10 m at 90 percent confidence level, both attributed to PRL 136, 120201 (2026). The two statements differ in both the value and the confidence level. Nothing in this paper's argument depends on which is correct, because the harness verdicts are internal to the models tested, but any future use of this bound in a quantitative claim must first resolve which statement the source supports.) The simple ladder lifetimes in the same record (Section 4.6) indicate the dynamic range involved: the Diósi-Penrose structure assigns lifetimes from 7.3529e+23 s for an electron superposition down to 9.7247e-27 s for a cat-scale superposition, which is why the passing Tilloy-Diosi form of Section 4.4 is a statement about a very specific, already-named theoretical structure and not about a generic new coupling.

The relevance of this context is directional. Models that fail K1 to K3 fail before experiment is consulted, at zero marginal cost. Models that pass do so by being known theory, whose experimental status is already settled: continuous measurement, confirmed wherever monitored systems are described, and the collapse family, constrained by bounds such as those above. No candidate of the forms tested here is both novel and unconstrained, and that is what the sieve is for: it tells the researcher which of those two words, novel or viable, has to be given up.

## 7 Discussion

The harness exists to be used before interpretation, by anyone holding a newly written coupling term. The workflow is short. Write the candidate down exactly. Grade it as written against K1, K2 and K3 at fixed settings. If it fails, the failure mode is itself the diagnosis: frequencies pinned at 0 and 1 mean the coupling does not move populations in its own basis; frequencies displaced uniformly mean a drift has broken the martingale; frequencies compressed toward 1/2 mean the coupling averages away before selection completes; a large K3 failure on the tested decomposition indicates that the localisation was bought with ensemble nonlinearity in that test. If it passes, the next step is not celebration but identification: compare the passing equations, term by term, against the standard theories they might secretly be. Section 5 is the worked example, and its outcome should be expected often enough that checking for it belongs in the protocol.

The limits of the harness should be stated as plainly as its results. It covers one qubit; multipartite systems, where K3's signalling content becomes operational rather than diagnostic, are not implemented. It covers Markovian diffusive dynamics; jump unravellings and non-Markovian records are not implemented, and a candidate whose natural form is a jump process would need its K1 analogue stated with the corresponding innovation structure, a generalisation we expect to be routine but have not carried out. The K2 SI conversion assumes a 1 micro-eV splitting to fix units, and the CSL figures assume rC = 1e-7 m and the GRW and Adler reference values; both are stated inputs, not findings. The underlying harness report declines to reduce the full lambda against rC bound curve to one number; readers who need the curve should go to Bassi et al. (2013) and subsequent bound compilations. Finally, the identity result is an identity for one qubit and one monitored observable. It does not prove that every conceivable observational coupling is continuous measurement in disguise; it proves that this one is, and it supplies the method, term substitution plus harness grading, by which any other can be checked.

The verdict language has a corresponding scope. A K1 pass for a stochastic model means that no departure from the Born frequencies was detected at the five tested populations with the stated trajectory counts and intervals; it does not prove exact Born invariance. Where the Born result follows from a martingale property of the model's own population equation, the run checks that the code implements that construction. A K3 pass for Model A or a candidate diffusion likewise concerns the specified decomposition and sampling floor, not all density operators, decompositions, or entangled extensions. Model B is an exception in the opposite direction: its ensemble linearity is exact because its dynamics is unitary conjugation. K2 is narrower still: it audits unconditional mean energy under the stated Hamiltonian. Heating can reject a model that promises isolated-system conservation or a contradicted energetic role, but heating alone does not reject an effective open-system model, and acceptable mean energy behaviour does not prove a complete conservation ledger.

We also record what this study did not find, because it is the finding most likely to be misreported. No candidate tested here produced a discrepancy with standard quantum mechanics: every failure was a failure to reproduce it, and every pass was a reproduction of it. A falsification harness that only ever falsifies the candidates is doing its job, and its job is worth doing in the open, with frozen code and rerunnable numbers, so that the next candidate can be graded in an afternoon rather than debated for a season.

## 8 Reproducibility statement

All numerical results in this paper come from frozen runs whose code, seeds, and machine-readable outputs are preserved in the research record of this study. The kill-test harness uses fixed seeds (20261118 to 20261122 for the original runs; 20270110 to 20270119 for the identity run), the numpy PCG64 generator, Euler-Maruyama integration at dt = 0.005 for diffusive trajectories, and exact diagonalisation for deterministic candidates. Every verdict reported here was reproduced in an independent regrade executed from byte-identical copies of the code on a different machine, as detailed in Sections 3.3 and 4.6: simulated frequencies matched exactly, and all other differences were at or below 3.1e-16 in analytic floating-point quantities. A second execution of the identity run's runner reproduced its results file with no differing byte. Section 3.3 gives a protocol by which an outside group can reproduce the verdicts from this paper alone, without access to the original code. Where the record contains two statements of the same quantity (the Model C trace distance at two precisions, and the XENONnT bound in two forms), both statements are reported in the text rather than resolved silently.

## Appendix: The Identity, Stated Precisely

This appendix separates the algebraic identity from the numerical consistency checks. The numerical runs show that the implementation behaves as the equations require at the stated sampling scales. They do not prove the identity. The identity is the term substitution in Section A.3, under the assumptions in Section A.4, and its level is fixed in Section A.5.

### A.1 The failed coupling and the repaired form

The failed coupling is Model B,

V = -W |1><1|, with W = 1 constant,

added to the Hamiltonian and evolved unitarily. In the sigma_z basis the term is diagonal. It changes a phase and leaves every population unchanged, so its population equation is dp = 0. It contains no stochastic term and no measurement record. Those facts, not a sampling fluctuation, are the source of its exact K1 failure.

The repaired ansatz is the conditional stochastic master equation

d rho_c = -i [H, rho_c] dt + gamma D[sigma_z] rho_c dt + sqrt(eta gamma) H[sigma_z] rho_c dW_t,

where

D[A] rho = A rho A^dagger - (1/2) {A^dagger A, rho}

and

H[A] rho = A rho + rho A^dagger - Tr[(A + A^dagger) rho] rho.

The stochastic equations in this appendix are in Ito form, with dW_t a Wiener increment. Averaging the conditional equation over the noise gives the unconditional equation

d rho/dt = -i [H, rho] + gamma D[sigma_z] rho.

At eta = 1, the corresponding pure-state equation is

d|psi_c> = [-iH - (gamma/2)(sigma_z - <sigma_z>)^2] |psi_c> dt + sqrt(gamma) (sigma_z - <sigma_z>) |psi_c> dW_t.

In the special setting [H, sigma_z] = 0, the conditional population obeys

dp = -4 sqrt(eta gamma) p (1 - p) dW_t.

The source derivation also records the ensemble heating functional

dE/dt = gamma Tr[(sigma_z H sigma_z - H) rho].

The heating functional is a consequence of the ensemble equation in the stated setting. It is not an additional term in the identity map.

The measurement record requires a qualification. The identity record underlying this paper identifies the record dY in the repair model as a homodyne-type diffusive record whose conditional mean is set by <c + c^dagger>_c, and it identifies the overall record gain as a labelling choice. It does not reproduce the record equation itself. This appendix therefore does not quote a symbol-by-symbol record equation. That is an explicit limitation: the state-equation identity is displayed in full, while the record is identified by its class, its conditional mean, and its gain status, not by a reproduced formula.

### A.2 Reference equation

The reference form used for the map is the diffusive stochastic master equation in Wiseman and Milburn, Quantum Measurement and Control, for a monitored channel with Lindblad operator c and efficiency eta:

d rho_c = -i [H, rho_c] dt + D[c] rho_c dt + sqrt(eta) H[c] rho_c dW_t.

The same H appears in the reference equation and in the ansatz. The map does not replace H, add a term to it, or remove a term from it.

### A.3 Parameter substitution and numerical factors

Set

c := sqrt(gamma) sigma_z.

For real gamma >= 0,

D[sqrt(gamma) A] rho = gamma D[A] rho.

For a nonnegative real factor and Hermitian A, the innovation superoperator is homogeneous of degree one in its operator argument:

H[sqrt(gamma) A] rho = sqrt(gamma) H[A] rho.

Applying the two identities with A = sigma_z gives

D[c] rho_c = gamma D[sigma_z] rho_c

and

sqrt(eta) H[c] rho_c dW_t = sqrt(eta gamma) H[sigma_z] rho_c dW_t.

Substitution therefore converts the reference equation into the ansatz equation symbol for symbol. The Hamiltonian drift, the dissipator, and the innovation term map individually. Averaging after substitution gives the unconditional Lindblad equation displayed in Section A.1. At eta = 1, the same substitution gives the pure-state stochastic Schroedinger equation displayed there.

The factors have distinct roles and are not interchangeable conventions. The coefficient gamma in gamma D[sigma_z] corresponds to a sigma_z coherence dephasing rate of 2 gamma in the ensemble equation. The efficiency eta is the observed fraction of the one channel: the standard split writes the channel in observed and unobserved parts proportional to sqrt(eta) c and sqrt(1 - eta) c, respectively. In the population equation, the source derivation records the equivalent convention gamma_J = 4 gamma, under which

dp = -4 sqrt(eta gamma) p (1 - p) dW_t = +2 sqrt(eta gamma_J) p (1 - p) dW'_t,

with the sign absorbed into the Wiener-process convention. A Stratonovich presentation is a presentation of the same equation only after the corresponding drift conversion; it is not a separate unconverted identity claim.

### A.4 Assumptions

The identity in Section A.5 is stated under all of the following assumptions.

1. The system is one qubit, with state space C^2, sigma_z = diag(+1, -1), and hbar = 1.
2. There is one monitored diffusive channel, with Hermitian channel operator proportional to sigma_z.
3. The same Hamiltonian H and the same initial conditional state are used on both sides of the substitution. The source derivation does not select a special initial state for the stochastic master equation identity.
4. The strength satisfies gamma >= 0, and the efficiency satisfies 0 <= eta <= 1.
5. The noise is a Markovian white Wiener process, and the stochastic equations are read in the Ito sense stated in Section A.1.
6. The stochastic master equation applies to conditional density operators, including mixed conditional states. The displayed pure-state stochastic Schroedinger equation is the eta = 1 case and presupposes a pure conditional state.
7. The SME identity does not require [H, sigma_z] = 0. The driftless population equation in Section A.1 is a consequence in the additional commuting setting.
8. The ansatz does not derive the monitored sigma_z channel from a more microscopic system and apparatus Hamiltonian. The source identity record marks that absence as an open item. The identity therefore begins after the monitored channel has been specified.

### A.5 The precise identity claimed

The identity claimed is equality of stochastic differential equations after parameter substitution.

More precisely, after c := sqrt(gamma) sigma_z, and under the assumptions in Section A.4, the repaired ansatz has exactly the Wiseman and Milburn diffusive stochastic master equation for the monitored channel. At eta = 1, it also has exactly the corresponding pure-state stochastic Schroedinger equation. Equality of the unconditional master equations follows by averaging and is a consequence of the SDE identity, not the level at which the claim stops. Equality in distribution of conditional trajectories and equality of terminal outcome probabilities also follow when the same initial condition and Wiener law are imposed, but neither is the primary claim.

In quantum probability terminology, the source identity record states that this diffusive equation, for one Hermitian observed channel, is the Belavkin filtering equation. No numerical value in this paper is taken from a Belavkin source, and the displayed derivation is the Wiseman and Milburn presentation.

The claim does not extend by itself to jump unravellings, non-Markovian records, other monitored observables, higher-dimensional systems, or a microscopic derivation of the channel. It is not an identity with CSL: the structural comparison in Section 5.2 remains a comparison, not a map. It is also not a minimality claim. No admissible class of repairs and no criterion for choosing a least repair are defined in the sources, so the text claims only that the repair displayed here satisfies the stated requirements and is identical, at the level fixed in this section, to standard continuous measurement.

## References

Bassi, A., Lochan, K., Satin, S., Singh, T. P., and Ulbricht, H. (2013). Models of wave-function collapse, underlying theories, and experimental tests. Reviews of Modern Physics, 85, 471. DOI: 10.1103/RevModPhys.85.471.

Donadi, S., et al. (2021). Nature Physics, 17, 74-78. DOI: 10.1038/s41567-020-1008-4.

Ghirardi, G. C., Pearle, P., and Rimini, A. (1990). Markov processes in Hilbert space and continuous spontaneous localization of systems of identical particles. Physical Review A, 42, 78. DOI: 10.1103/PhysRevA.42.78.

Gisin, N. (1989). Stochastic quantum dynamics and relativity. Helvetica Physica Acta, 62, 363.

Wiseman, H. M., and Milburn, G. J. (2010). Quantum Measurement and Control. Cambridge University Press. DOI: 10.1017/CBO9780511813948.

Belavkin, V. P. (1992). Communications in Mathematical Physics, 146, 611-635. DOI: 10.1007/BF02097018.

Tilloy, A., and Diosi, L. (2016). Sourcing semiclassical gravity from spontaneously localized quantum matter. Physical Review D, 93, 024026. DOI: 10.1103/PhysRevD.93.024026. (Bibliographic details confirmed against Crossref.)

XENONnT Collaboration (2026). Physical Review Letters, 136, 120201. See the footnote in Section 6 on the two preserved statements of this bound.

## AI assistance disclosure

AI tools were used in the preparation of this work: they assisted in drafting and revising the text of this manuscript, in reviewing the simulation code and run records on which the reported numbers rest, and in checking the internal consistency of the reported values against the frozen results files. All scientific decisions, including the design of the tests, the grading criteria, and the decision to report the preserved discrepancies noted in Sections 4.3 and 6, are the author's. The author has reviewed the final text and takes full responsibility for its contents, including every number quoted. AI tools are not authors of this work and are not listed as such.
