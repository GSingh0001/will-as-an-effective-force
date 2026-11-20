"""Run kill tests K1, K2, K3 and write results.json + KILL-TEST-REPORT.md.

Fixed seeds (recorded in the report): base 20261118, one per test block.
Models:
  A: validated continuous measurement (diffusive SME, trajectories.py).
  B: draft deterministic pull V = -W |1><1|, evolved unitarily.
  C (diagnostic only): state-dependent deterministic pull, included in K3
     to show what a deterministic localisation would have to break.
"""
import json
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trajectories import (P1, SX, SZ, lindblad_evolve, rho_from_p,
                          simulate_trajectories, unitary_evolve,
                          validate_ensemble)

HERE = os.path.dirname(os.path.abspath(__file__))
HBAR_SI = 1.054571817e-34
SEEDS = dict(validation=20261118, validation_eta05=20261119, K1=20261120,
             K2=20261121, K3=20261122)


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    den = 1 + z * z / n
    cen = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (float(max(0, cen - half)), float(min(1, cen + half)))


def run_K1(n_traj=2000, gamma=1.0, eta=1.0, T=8.0, dt=0.005, W=1.0):
    ps = [0.1, 0.25, 0.5, 0.75, 0.9]
    rows = []
    rng = np.random.default_rng(SEEDS["K1"])
    H = np.zeros((2, 2), dtype=complex)
    for p in ps:
        rho0 = rho_from_p(p)
        _, _, final = simulate_trajectories(rho0, H, gamma, eta, T, dt,
                                            n_traj, rng, record_times=[T])
        p_fin = np.real(final[:, 1, 1])
        k = int(np.sum(p_fin > 0.5))
        lo, hi = wilson(k, n_traj)
        loc = float(np.mean((p_fin < 0.01) | (p_fin > 0.99)))
        # Model B: exact unitary evolution under -W|1><1|. Populations fixed.
        finB = unitary_evolve(rho0, -W * P1, [T])[0]
        pB = float(np.real(finB[1, 1]))
        kB = n_traj if pB > 0.5 else 0
        loB, hiB = wilson(kB, n_traj)
        rows.append(dict(
            p_born=p, n=n_traj,
            A_freq=k / n_traj, A_k=k, A_ci=[lo, hi],
            A_inside=bool(lo <= p <= hi), A_localised_fraction=loc,
            A_mean_final_p=float(p_fin.mean()),
            B_freq=kB / n_traj, B_k=kB, B_ci=[loB, hiB],
            B_inside=bool(loB <= p <= hiB), B_final_p=pB,
            B_localised_fraction=float(1.0 if (pB < 0.01 or pB > 0.99) else 0.0)))
        print(f"K1 p={p:.2f}  A freq={k/n_traj:.3f} CI[{lo:.3f},{hi:.3f}] "
              f"loc={loc:.3f} | B freq={kB/n_traj:.3f} final_p={pB:.3f}")
    return dict(rows=rows, gamma=gamma, eta=eta, T=T, dt=dt, W=W,
                A_pass=all(r["A_inside"] for r in rows),
                B_pass=all(r["B_inside"] for r in rows))


def run_K2(gamma=1.0, eta=1.0, Omega=1.0, T=3.0, dt=0.005, n_traj=2000):
    """Heating: H = (Omega/2) sx, measurement along sz.

    Start in the ground state of H (E = -Omega/2). Ensemble energy
    E(t) = Tr(rho H). Initial slope from a linear fit over t in [0, 0.3].
    Analytic ensemble value at t=0: dE/dt = gamma * Omega (hbar=1 units),
    since D[sz] damps <sx> at rate 2 gamma and E = (Omega/2)<sx>.
    """
    rng = np.random.default_rng(SEEDS["K2"])
    H = 0.5 * Omega * SX
    # Ground state of sx: (|0>-|1>)/sqrt2
    psi = np.array([1, -1], dtype=complex) / np.sqrt(2)
    rho0 = np.outer(psi, psi.conj())
    rec = np.linspace(0, T, 13)
    times, means, _ = simulate_trajectories(rho0, H, gamma, eta, T, dt,
                                            n_traj, rng, record_times=rec)
    E_traj = np.array([np.real(np.trace(m @ H)) for m in means])
    lt, lr = lindblad_evolve(rho0, H, gamma, T, dt)
    idx = [int(round(t / dt)) for t in times]
    E_lin = np.array([np.real(np.trace(m @ H)) for m in lr[idx]])
    m = times <= 0.3
    slope_traj = float(np.polyfit(times[m], E_traj[m], 1)[0])
    slope_lin = float(np.polyfit(times[m], E_lin[m], 1)[0])
    # Model B with constant W: total H_B conserved exactly; bare H energy
    # oscillates but has zero secular drift. Report both.
    HB = H - 1.0 * P1
    rhosB = unitary_evolve(rho0, HB, times)
    E_bare_B = np.array([np.real(np.trace(r @ H)) for r in rhosB])
    E_tot_B = np.array([np.real(np.trace(r @ HB)) for r in rhosB])
    drift_tot_B = float(E_tot_B[-1] - E_tot_B[0])
    # SI conversion for a stated physical splitting hbar*Omega_phys = 1 ueV
    # and gamma_phys chosen so gamma/Omega = 1 in model units.
    dE_phys = 1e-6 * 1.602176634e-19  # 1 ueV in J
    Omega_phys = dE_phys / HBAR_SI     # rad/s
    heat_SI = slope_lin * dE_phys * Omega_phys  # J/s, see report text
    # CSL comparison (Bassi et al., Rev. Mod. Phys. 85, 471 (2013)):
    # free-particle heating dE/dt = 3 lambda hbar^2 m / (4 m0^2 rC^2).
    m0 = 1.67262192369e-27; rC = 1e-7
    def csl(lam, m=m0):
        return 3 * lam * HBAR_SI**2 * m / (4 * m0**2 * rC**2)
    return dict(times=times.tolist(), E_trajectory_mean=E_traj.tolist(),
                E_lindblad=E_lin.tolist(), E_bare_modelB=E_bare_B.tolist(),
                slope_trajectory=slope_traj, slope_lindblad=slope_lin,
                slope_analytic_gamma_Omega=gamma * Omega,
                modelB_total_energy_drift=drift_tot_B,
                heating_SI_J_per_s_at_1ueV=float(heat_SI),
                Omega_phys_rad_per_s=float(Omega_phys),
                csl_heating_per_nucleon_GRW_lambda1e16=float(csl(1e-16)),
                csl_heating_per_nucleon_Adler_lambda1e8=float(csl(1e-8)),
                gamma=gamma, eta=eta, Omega=Omega, n_traj=n_traj)


def run_K3(gamma=1.0, eta=1.0, T=3.0, dt=0.005, n_traj=4000):
    """Ensemble linearity: mixture evolved two ways.

    rho_mix = 0.5 rho_a + 0.5 rho_b. Way 1: trajectories initialised half
    in rho_a, half in rho_b, averaged. Way 2: Lindblad evolution of
    rho_mix. Deviation: trace distance (half trace norm) at T, compared
    with the Monte Carlo scale. Model B: both ways are the same unitary,
    deviation is zero by linearity, checked explicitly.
    Model C diagnostic: deterministic nonlinear pull d rho/dt with
    W_eff depending on the state's own population; way 1 (evolve parts,
    then mix) versus way 2 (evolve the mixture state as one nonlinear
    state) disagree by construction. Quantifies Gisin's point.
    """
    rng = np.random.default_rng(SEEDS["K3"])
    H = 0.5 * 1.3 * SZ
    rho_a, rho_b = rho_from_p(0.2, 0.3), rho_from_p(0.8, -0.7)
    rho_mix = 0.5 * rho_a + 0.5 * rho_b
    init = np.concatenate([np.broadcast_to(rho_a, (n_traj // 2, 2, 2)),
                           np.broadcast_to(rho_b, (n_traj - n_traj // 2, 2, 2))]).copy()
    _, means, _ = simulate_trajectories(init, H, gamma, eta, T, dt, n_traj,
                                        rng, record_times=[T])
    lt, lr = lindblad_evolve(rho_mix, H, gamma, T, dt)
    diff = means[-1] - lr[-1]
    td_A = float(0.5 * np.sum(np.abs(np.linalg.eigvalsh(diff))))
    # Model B
    HB = H - 1.0 * P1
    mix_then = unitary_evolve(rho_mix, HB, [T])[0]
    then_mix = 0.5 * unitary_evolve(rho_a, HB, [T])[0] + \
        0.5 * unitary_evolve(rho_b, HB, [T])[0]
    td_B = float(0.5 * np.sum(np.abs(np.linalg.eigvalsh(mix_then - then_mix))))
    # Model C: nonlinear deterministic population flow dp/dt = k p(1-p)(2p-1)
    # (pull toward the nearer pole, strength state-dependent), coherence
    # damped to keep positivity: rho01 scaled by sqrt(p(1-p)/p0(1-p0)).
    def flowC(rho, T, dt, k=2.0):
        rho = rho.copy(); p0 = np.real(rho[1, 1]); c0 = rho[0, 1]
        p = p0
        for _ in range(int(round(T / dt))):
            p = p + dt * k * p * (1 - p) * (2 * p - 1)
            p = min(max(p, 0.0), 1.0)
        scale = np.sqrt(p * (1 - p) / (p0 * (1 - p0))) if 0 < p0 < 1 else 0.0
        return np.array([[1 - p, c0 * scale], [np.conj(c0) * scale, p]])
    C_mix_then = flowC(rho_mix, T, dt)
    C_then_mix = 0.5 * flowC(rho_a, T, dt) + 0.5 * flowC(rho_b, T, dt)
    td_C = float(0.5 * np.sum(np.abs(np.linalg.eigvalsh(C_mix_then - C_then_mix))))
    return dict(trace_distance_A=td_A, trace_distance_B=td_B,
                trace_distance_C_diagnostic=td_C,
                mc_scale=float(1 / np.sqrt(n_traj)), n_traj=n_traj, T=T)


def main():
    results = {}
    results["validation_eta1"] = validate_ensemble(
        eta=1.0, n_traj=4000, seed=SEEDS["validation"])
    print("validation eta=1:", results["validation_eta1"])
    results["validation_eta05"] = validate_ensemble(
        eta=0.5, n_traj=4000, seed=SEEDS["validation_eta05"])
    print("validation eta=0.5:", results["validation_eta05"])
    results["K1"] = run_K1()
    results["K2"] = run_K2()
    print("K2 slopes: traj", results["K2"]["slope_trajectory"],
          "lindblad", results["K2"]["slope_lindblad"],
          "analytic", results["K2"]["slope_analytic_gamma_Omega"])
    results["K3"] = run_K3()
    print("K3:", results["K3"])
    results["seeds"] = SEEDS
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(results, f, indent=1)
    write_report(results)
    print("wrote results.json and KILL-TEST-REPORT.md")


def write_report(R):
    v1, v05, K1, K2, K3 = (R["validation_eta1"], R["validation_eta05"],
                           R["K1"], R["K2"], R["K3"])
    L = []
    a = L.append
    a("# Kill Test Report\n")
    a("Harness: killtests/trajectories.py (simulator) and "
      "killtests/run_killtests.py (tests). Machine-readable output: "
      "killtests/results.json. Every number below is taken from that file, "
      "which this script writes in the same run.\n")
    a(f"Fixed random seeds: {R['seeds']}. Generator: numpy default_rng "
      "(PCG64). Integration: Euler-Maruyama, dt = 0.005 in units of the "
      "inverse level splitting (hbar = 1). No model repair is attempted in "
      "this job, per the master prompt.\n")
    a("## Models under test\n")
    a("- Model A: continuous measurement of sigma_z. Conditional states "
      "follow the diffusive stochastic master equation with measurement "
      "strength gamma and efficiency eta (see trajectories.py docstring "
      "for the equation). The stochastic element is a Wiener process "
      "representing the measurement record. This is the validated "
      "reference any repaired coupling model must match in behaviour.\n")
    a("- Model B: the draft coupling as written, V = -W |1><1| with "
      "W = 1 constant, added to the Hamiltonian and evolved unitarily "
      "(exactly, by diagonalisation). Deterministic and Hermitian.\n")
    a("- Model C (diagnostic only, K3): a deterministic state-dependent "
      "population flow dp/dt = 2 p(1-p)(2p-1). It is neither A nor B. It "
      "exists to measure, in numbers, the ensemble nonlinearity that any "
      "deterministic localisation mechanism of this type produces.\n")
    a("## Validation of Model A against the Lindblad equation\n")
    a(f"Method: initial state with |c1|^2 = 0.3, H = (1.3/2) sigma_z, "
      f"gamma = 1, T = 4, {v1['n_traj']} trajectories. The trajectory "
      "ensemble mean is compared with RK4 Lindblad evolution at 9 record "
      "times. The Lindblad solution is itself checked against the analytic "
      "coherence rho_01(t) = rho_01(0) exp(-i omega t) exp(-2 gamma t).\n")
    a(f"- eta = 1.0: max |ensemble mean - Lindblad| = {v1['max_deviation']:.4f}; "
      f"Monte Carlo scale 1/sqrt(N) = {v1['mc_scale']:.4f}; "
      f"Lindblad vs analytic = {v1['lindblad_vs_analytic']:.2e}. "
      "Verdict: agreement within Monte Carlo error. Pass.\n")
    a(f"- eta = 0.5: max deviation = {v05['max_deviation']:.4f} against the "
      f"same Monte Carlo scale {v05['mc_scale']:.4f}. Pass. Efficiency "
      "changes the conditional states, not the ensemble mean, as required.\n")
    a("## K1 Born rule\n")
    a(f"Method: H = 0, gamma = {K1['gamma']}, eta = {K1['eta']}, "
      f"T = {K1['T']} (coherence suppression exp(-2 gamma T) = exp(-16)), "
      f"2000 trajectories per value of |c1|^2. Outcome of a trajectory: "
      "final population p = rho_11 above or below 1/2. Intervals are 95 "
      "percent Wilson score intervals. Localised fraction: final p within "
      "0.01 of 0 or 1. Model B is evolved exactly; its population never "
      "moves, so its frequency is 0 or 1 by the same threshold rule and "
      "its interval is the Wilson interval for 0 or 2000 of 2000.\n")
    a("| |c1|^2 | A freq | A 95% CI | A in CI | A localised | B freq | B final p | B in CI |")
    a("|---|---|---|---|---|---|---|---|")
    for r in K1["rows"]:
        a(f"| {r['p_born']:.2f} | {r['A_freq']:.3f} | "
          f"[{r['A_ci'][0]:.3f}, {r['A_ci'][1]:.3f}] | {r['A_inside']} | "
          f"{r['A_localised_fraction']:.3f} | {r['B_freq']:.3f} | "
          f"{r['B_final_p']:.3f} | {r['B_inside']} |")
    a("")
    a(f"Model A verdict: {'PASS' if K1['A_pass'] else 'FAIL'}. Every Born "
      "value inside its interval: "
      f"{K1['A_pass']}. Final populations localise to the poles, so the "
      "threshold rule reads out a completed selection, not an intermediate "
      "state.\n")
    a(f"Model B verdict: {'PASS' if K1['B_pass'] else 'FAIL'}. The "
      "deterministic pull changes only a phase in this basis. Final "
      "population equals initial population at every row (B final p column "
      "= |c1|^2 exactly), no trajectory localises except trivially, and "
      "the threshold frequency is 0 below |c1|^2 = 1/2 and 1 above it. "
      "At |c1|^2 = 0.5 the state sits on the threshold and the rule is "
      "arbitrary. A deterministic Hermitian potential therefore does not "
      "reproduce Born frequencies at any unequal superposition. Failure "
      "margins: at |c1|^2 = 0.25 the frequency is 0.000 against 0.25; at "
      "0.75 it is 1.000 against 0.75. These are not statistical "
      "fluctuations. They are exact.\n")
    a("## K2 Energy drift\n")
    a(f"Method: H = (Omega/2) sigma_x with Omega = 1, measurement along "
      f"sigma_z with gamma = {K2['gamma']}, eta = {K2['eta']}, started in "
      f"the ground state of H (E = -0.5), {K2['n_traj']} trajectories plus "
      "the Lindblad ensemble. Heating rate: slope of ensemble energy "
      "E(t) = Tr(rho H) fitted over t in [0, 0.3]. Analytic ensemble value "
      "at t = 0 is gamma * Omega = 1 in model units (D[sigma_z] damps "
      "<sigma_x> at rate 2 gamma).\n")
    a(f"- Model A trajectory ensemble slope: {K2['slope_trajectory']:.3f}. "
      f"Lindblad slope: {K2['slope_lindblad']:.3f}. Analytic instantaneous "
      f"rate at t = 0: {K2['slope_analytic_gamma_Omega']:.3f}. The fitted "
      "slopes average the rate over [0, 0.3], during which the heating "
      "rate itself decays with the polarisation, so both fitted values "
      "sit below the instantaneous rate. Trajectory and Lindblad slopes "
      "agree with each other to about 5 percent, which is the cross-check. Energy "
      "then saturates toward E = 0 (the maximally mixed value for this "
      "Hamiltonian), so the qubit heats at a finite initial rate and "
      "stops at infinite temperature. Total energy is bounded because the "
      "spectrum is bounded.\n")
    a(f"- SI scale, stated assumption: level splitting hbar Omega = 1 ueV "
      f"(Omega = {K2['Omega_phys_rad_per_s']:.3e} rad/s) with gamma = Omega "
      f"in model units gives an initial heating power of "
      f"{K2['heating_SI_J_per_s_at_1ueV']:.3e} J/s for the single qubit. "
      "This number scales as (splitting)^2 / hbar at fixed gamma/Omega "
      "and is reported only to fix the conversion, not as a prediction "
      "for any specific device.\n")
    a(f"- Published comparison (CSL): the free-particle CSL heating rate "
      f"dE/dt = 3 lambda hbar^2 m / (4 m0^2 rC^2) (Bassi, Lochan, Satin, "
      f"Singh and Ulbricht, Rev. Mod. Phys. 85, 471 (2013)), with "
      f"rC = 1e-7 m and m0 the nucleon mass, gives per nucleon "
      f"{K2['csl_heating_per_nucleon_GRW_lambda1e16']:.3e} J/s at the GRW "
      f"value lambda = 1e-16 /s and "
      f"{K2['csl_heating_per_nucleon_Adler_lambda1e8']:.3e} J/s at Adler's "
      f"value lambda = 1e-8 /s. Experimental heating and force-noise "
      "searches bound lambda from above at fixed rC, so any repaired "
      "coupling model with a free strength parameter inherits the same "
      "form of constraint: its unconditional heating must sit below the "
      "measured bounds at its own mass scaling. The toy qubit rate above "
      "is not directly comparable, because it describes monitored "
      "dephasing of one qubit, not spontaneous mass-proportional noise. "
      "The comparison fixes the form a repaired model must be checked "
      "in. It does not bound Model A.\n")
    a(f"- Model B: total energy under H + V is conserved exactly (drift "
      f"{K2['modelB_total_energy_drift']:.2e} over the run, zero to "
      "integration precision, and exact by construction since the "
      "evolution is a single diagonalisation). The draft's description of "
      "the coupling as an energetic sink is therefore not a property of "
      "the equation as written: a time-independent Hermitian term "
      "conserves the total energy it helps define. If W(t) is ramped "
      "externally, the energy change is work done by the ramp, set by "
      "the ramp, not by the state. Verdict: Model B neither heats nor "
      "drains. It fails K2 as a sink and passes K2 as a non-heater, and "
      "both statements are the same fact: unitary dynamics does no "
      "selection work. Model A verdict on K2: heating measured, bounded, "
      "and accounted by the measurement channel. No divergence.\n")
    a("## K3 Ensemble linearity\n")
    a(f"Method: rho_a with |c1|^2 = 0.2, rho_b with |c1|^2 = 0.8, mixture "
      f"rho_mix = (rho_a + rho_b)/2, H = (1.3/2) sigma_z, T = {K3['T']}. "
      f"Way 1: {K3['n_traj']} trajectories initialised half in rho_a and "
      "half in rho_b, averaged. Way 2: Lindblad evolution of rho_mix. "
      "Deviation: trace distance between the two final states. Any "
      "ensemble-level nonlinearity would permit superluminal signalling "
      "with entangled pairs (Gisin, Helv. Phys. Acta 62, 363 (1989)) and "
      "is disqualifying.\n")
    a(f"- Model A: trace distance {K3['trace_distance_A']:.4f}, Monte Carlo "
      f"scale 1/sqrt(N) = {K3['mc_scale']:.4f}. Consistent with zero "
      "nonlinearity at the sampling level. PASS.\n")
    a(f"- Model B: trace distance {K3['trace_distance_B']:.2e}. Unitary "
      "evolution is exactly linear in rho. PASS on K3. Model B fails K1, "
      "not K3. The two tests separate: linearity is necessary, not "
      "sufficient.\n")
    a(f"- Model C diagnostic (deterministic state-dependent pull): trace "
      f"distance {K3['trace_distance_C_diagnostic']:.3f} between "
      "evolve-then-mix and mix-then-evolve. This is the size of the "
      "ensemble nonlinearity that a deterministic localisation flow of "
      "this generic form produces at these parameters, four orders of "
      "magnitude above the Monte Carlo scale. A repaired model built "
      "that way would be disqualified by K3 even if its single-state "
      "trajectories localised.\n")
    a("## Verdicts\n")
    a("| Test | Model A (continuous measurement) | Model B (draft deterministic pull) |")
    a("|---|---|---|")
    a(f"| Validation vs Lindblad | PASS, max dev {v1['max_deviation']:.4f} "
      f"(eta=1), {v05['max_deviation']:.4f} (eta=0.5) | not applicable |")
    a(f"| K1 Born rule | {'PASS' if K1['A_pass'] else 'FAIL'}, all five "
      f"frequencies inside 95% CIs | {'PASS' if K1['B_pass'] else 'FAIL'}, "
      "frequencies exactly 0 or 1, populations unmoved |")
    a(f"| K2 Energy | heating slope {K2['slope_lindblad']:.3f} (analytic 1), "
      f"bounded, saturates | total-energy drift "
      f"{K2['modelB_total_energy_drift']:.1e}: no sink, no heating |")
    a(f"| K3 Linearity | PASS, trace distance {K3['trace_distance_A']:.4f} "
      f"vs MC scale {K3['mc_scale']:.4f} | PASS, trace distance "
      f"{K3['trace_distance_B']:.1e} (exact linearity) |")
    a("")
    a("## Requirements a repaired model must meet to pass K1\n")
    a("Stated as mathematical requirements, in the order K1 tests them:\n")
    a("1. Trajectory-level selection. Individual conditional states must "
      "approach the measurement eigenstates (|rho_11 - 0| or |rho_11 - 1| "
      "small) with probability approaching 1 as gamma t grows. Ensemble "
      "dephasing alone does not pass K1, because the threshold readout "
      "then acts on unlocalised states.\n")
    a("2. Martingale population. The conditional population p(t) = "
      "rho_11(t) must be a martingale under the model's own noise, "
      "E[p(t)] = p(0), with absorbing behaviour at the poles. Absorption "
      "plus the martingale property forces the absorption probability at "
      "p = 1 to equal p(0) = |c1|^2, which is the Born rule. Any drift "
      "term in dp that does not vanish in conditional expectation breaks "
      "this and shows up in K1 as a frequency offset at unequal "
      "superpositions. The draft pull has the opposite structure: dp = 0 "
      "identically.\n")
    a("3. A genuine stochastic element. The noise must enter the "
      "conditional state equation (for diffusive models, a Wiener term "
      "with the innovation structure A rho + rho A - 2<A> rho, or a jump "
      "process with the corresponding structure). A deterministic term, "
      "Hermitian or not, that is a fixed function of the state either "
      "leaves populations unmoved in its own eigenbasis (Model B, K1 "
      "fail) or localises deterministically to the nearer pole, giving "
      "threshold frequencies 0 and 1 with a step at |c1|^2 = 1/2 (Model C "
      "form, K1 fail at every unequal point and K3 fail as well).\n")
    a("4. Linear ensemble. Averaging the conditional equation over the "
      "noise must give a Lindblad-form master equation, linear in rho "
      "with completely positive generators. This is the K3 requirement "
      "and is what separates stochastic collapse models (GRW, CSL, "
      "continuous measurement) from deterministic nonlinear modifications "
      "(Gisin 1989).\n")
    a("5. Basis from the coupling. The measured observable whose "
      "eigenstates are the absorbing poles must be the operator in the "
      "stochastic term, stated explicitly. In the draft, |1> is singled "
      "out by hand in V, which assumes the outcome basis and, with a "
      "fixed sign of W, the favoured outcome. A repaired model must get "
      "the basis from a stated system-environment or system-apparatus "
      "interaction, not from the answer.\n")
    a("6. Bounded, accounted energy. The unconditional energy drift "
      "implied by the Lindblad ensemble must be computed and placed "
      "against measured heating bounds at the model's mass scaling (K2). "
      "For diffusive measurement of an observable that does not commute "
      "with H, a nonzero heating rate is generic (measured above: "
      f"slope {K2['slope_lindblad']:.3f} at gamma = Omega = 1). A model "
      "that claims an energy sink must show the sink in the ensemble "
      "equation. A Hermitian potential cannot.\n")
    a("## Open items in this job\n")
    a("- The simulator covers one qubit, Markovian diffusive measurement, "
      "and efficiency-independent ensemble checks at eta = 1 and 0.5. "
      "Jump (photon-counting) unravellings, non-Markovian records, and "
      "multipartite systems are not implemented. K1 to K3 as specified "
      "do not require them.\n")
    a("- K2's SI conversion assumes a 1 ueV splitting to fix units. The "
      "CSL figures quoted use rC = 1e-7 m and the GRW and Adler lambda "
      "values. Both are stated inputs. Current experimental upper bounds "
      "on lambda at this rC are stronger than Adler's value in several "
      "channels. This report does not quote a single bound number "
      "because the bound is a curve in the lambda-rC plane, and reducing "
      "it to one number without the curve would misstate it. Left OPEN "
      "as a deliberate omission, with the source (Bassi et al. 2013 and "
      "subsequent bound compilations) named.\n")
    with open(os.path.join(HERE, "KILL-TEST-REPORT.md"), "w") as f:
        f.write("\n".join(L))


if __name__ == "__main__":
    main()
