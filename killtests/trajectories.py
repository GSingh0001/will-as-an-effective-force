"""Quantum trajectory simulator: two-level system under continuous measurement.

Model A (validated continuous measurement). Conditional state rho_c obeys
the diffusive stochastic master equation (Ito form, hbar = 1)

  d rho_c = -i [H, rho_c] dt + gamma D[sz] rho_c dt
            + sqrt(eta * gamma) H[sz] rho_c dW,

  D[A] rho = A rho A^dagger - (1/2) {A^dagger A, rho},
  H[A] rho = A rho + rho A^dagger - Tr[(A + A^dagger) rho] rho,

with A = sz (Hermitian, A^2 = I), measurement strength gamma >= 0 and
efficiency 0 <= eta <= 1. The ensemble average over the Wiener noise W
obeys the Lindblad master equation

  d rho / dt = -i [H, rho] + gamma D[sz] rho,

in which the coherence rho_01 decays as exp(-2 gamma t) (times the
Hamiltonian phase). validate_ensemble() checks this numerically.

Model B (draft deterministic pull). H_draft = H - W |1><1| with constant
W, evolved unitarily. Populations in the measurement basis are exactly
conserved when H is diagonal in that basis; no outcome is ever selected.

All simulations are vectorised over trajectories. Units: hbar = 1,
energies in units of the chosen level splitting, times in inverse
splitting. run_killtests.py converts the K2 heating rate to SI for a
stated physical splitting.
"""
from __future__ import annotations

import numpy as np

SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)
P1 = np.array([[0, 0], [0, 1]], dtype=complex)  # |1><1|


def rho_from_p(p: float, phase: float = 0.0) -> np.ndarray:
    """Pure state sqrt(1-p)|0> + exp(i phase) sqrt(p)|1> as density matrix."""
    psi = np.array([np.sqrt(1 - p), np.exp(1j * phase) * np.sqrt(p)], dtype=complex)
    return np.outer(psi, psi.conj())


def simulate_trajectories(rho0, H, gamma, eta, T, dt, n_traj, rng,
                          record_times=None):
    """Euler-Maruyama integration of the stochastic master equation.

    rho0 may be a single (2,2) state (broadcast to n_traj copies) or an
    (n_traj,2,2) array of per-trajectory initial states (used for K3
    mixtures). Returns (times, ensemble_mean[t], final_rhos) where times
    are the record times (default: every step).
    """
    rho0 = np.asarray(rho0, dtype=complex)
    if rho0.ndim == 2:
        rho = np.broadcast_to(rho0, (n_traj, 2, 2)).copy()
    else:
        rho = rho0.copy()
        n_traj = rho.shape[0]
    n_steps = int(round(T / dt))
    if record_times is None:
        record_steps = set(range(n_steps + 1))
    else:
        record_steps = {int(round(t / dt)) for t in record_times}
    times, means = [], []
    sq = np.sqrt(eta * gamma) if gamma > 0 else 0.0
    sqdt = np.sqrt(dt)
    for step in range(n_steps + 1):
        if step in record_steps:
            times.append(step * dt)
            means.append(rho.mean(axis=0))
        if step == n_steps:
            break
        # Deterministic part: -i[H,rho] + gamma (sz rho sz - rho)
        comm = H @ rho - rho @ H
        drho = -1j * comm + gamma * (SZ @ rho @ SZ - rho)
        # Stochastic part: sqrt(eta gamma) (sz rho + rho sz - 2<sz>rho) dW
        exp_sz = np.real(rho[:, 0, 0] - rho[:, 1, 1])
        innov = SZ @ rho + rho @ SZ - 2.0 * exp_sz[:, None, None] * rho
        dW = rng.standard_normal(n_traj) * sqdt
        rho = rho + drho * dt + sq * innov * dW[:, None, None]
        # Enforce Hermiticity and unit trace against integration drift.
        rho = 0.5 * (rho + rho.conj().transpose(0, 2, 1))
        tr = np.real(rho[:, 0, 0] + rho[:, 1, 1])
        rho = rho / tr[:, None, None]
    return np.array(times), np.array(means), rho


def lindblad_evolve(rho0, H, gamma, T, dt):
    """Deterministic Lindblad evolution by RK4. Returns (times, rhos)."""
    def rhs(rho):
        return -1j * (H @ rho - rho @ H) + gamma * (SZ @ rho @ SZ - rho)
    n = int(round(T / dt))
    rho = np.asarray(rho0, dtype=complex).copy()
    times, rhos = [0.0], [rho.copy()]
    for i in range(n):
        k1 = rhs(rho); k2 = rhs(rho + dt / 2 * k1)
        k3 = rhs(rho + dt / 2 * k2); k4 = rhs(rho + dt * k3)
        rho = rho + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        times.append((i + 1) * dt); rhos.append(rho.copy())
    return np.array(times), np.array(rhos)


def unitary_evolve(rho0, H_total, times):
    """Exact unitary evolution rho(t) = U rho0 U^dagger at given times."""
    w, V = np.linalg.eigh(H_total)
    out = []
    for t in times:
        U = V @ np.diag(np.exp(-1j * w * t)) @ V.conj().T
        out.append(U @ rho0 @ U.conj().T)
    return np.array(out)


def validate_ensemble(gamma=1.0, eta=1.0, omega=1.3, T=4.0, dt=0.005,
                      n_traj=4000, seed=20261118):
    """Compare trajectory ensemble mean with Lindblad evolution.

    Returns dict with max |mean - lindblad| over matrix elements and
    record times, plus the Monte Carlo scale 1/sqrt(n_traj) for scale.
    """
    rng = np.random.default_rng(seed)
    H = 0.5 * omega * SZ
    rho0 = rho_from_p(0.3, phase=0.4)
    rec = np.linspace(0, T, 9)
    times, means, _ = simulate_trajectories(rho0, H, gamma, eta, T, dt,
                                            n_traj, rng, record_times=rec)
    lt, lr = lindblad_evolve(rho0, H, gamma, T, dt)
    idx = [int(round(t / dt)) for t in times]
    dev = np.max(np.abs(means - lr[idx]))
    # Analytic coherence check: rho01(t) = rho01(0) exp(-i omega t) exp(-2 gamma t)
    ana01 = rho0[0, 1] * np.exp(-1j * omega * times) * np.exp(-2 * gamma * times)
    dev_ana = np.max(np.abs(lr[idx][:, 0, 1] - ana01))
    return dict(max_deviation=float(dev), lindblad_vs_analytic=float(dev_ana),
                mc_scale=float(1.0 / np.sqrt(n_traj)), n_traj=n_traj,
                gamma=gamma, eta=eta, dt=dt, T=T)
