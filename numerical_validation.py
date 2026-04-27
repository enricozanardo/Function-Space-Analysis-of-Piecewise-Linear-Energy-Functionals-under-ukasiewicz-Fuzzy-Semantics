"""Numerical validation of the analytical bounds derived in
'Function Space Analysis of Piecewise-Linear Energy Functionals
under Lukasiewicz Fuzzy Semantics'.

Reproducibility script accompanying Section 8 of the manuscript.

Usage
-----
    python numerical_validation.py [--output-dir results]

Dependencies: numpy >= 1.24, scipy >= 1.10, matplotlib >= 3.7.
Run-time: under 60 s on a single CPU core for the default settings.
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass, field

import numpy as np
from scipy.special import logsumexp

RNG = np.random.default_rng(20260426)


# ---------------------------------------------------------------------------
# Lukasiewicz energy and its gradient
# ---------------------------------------------------------------------------

@dataclass
class HornKB:
    """Horn-clause knowledge base: each clause is (head, body, weight).

    head :: int           index of the head atom
    body :: tuple[int]    indices of the body atoms (possibly empty)
    weight :: float       non-negative weight
    """
    n_atoms: int
    clauses: list[tuple[int, tuple[int, ...], float]] = field(default_factory=list)

    @property
    def W(self) -> float:
        return float(sum(abs(w) for _, _, w in self.clauses))

    # Lukasiewicz t-norm for a Horn clause body & implication
    @staticmethod
    def _luk_and(values: np.ndarray) -> np.ndarray:
        # max(0, sum(a_j) - (k-1))
        if values.size == 0:
            return np.ones(())
        return np.maximum(0.0, values.sum(axis=0) - (values.shape[0] - 1))

    def sigma(self, I: np.ndarray) -> np.ndarray:
        """Return the satisfaction vector sigma(I) of shape (m,) for I in [0,1]^n.

        I has shape (..., n_atoms); the leading axes are batch dimensions.
        """
        m = len(self.clauses)
        out = np.zeros(I.shape[:-1] + (m,))
        for k, (h, body, _) in enumerate(self.clauses):
            head = I[..., h]
            if body:
                body_vals = np.stack([I[..., j] for j in body], axis=0)
                body_truth = self._luk_and(body_vals)
            else:
                body_truth = np.ones_like(head)
            # Lukasiewicz implication: min(1, 1 - body + head)
            out[..., k] = np.minimum(1.0, 1.0 - body_truth + head)
        return out

    def energy(self, I: np.ndarray) -> np.ndarray:
        """E(I; w) = sum_i w_i sigma_i(I)."""
        s = self.sigma(I)
        w = np.array([w for _, _, w in self.clauses])
        return s @ w

    def grad_energy(self, I: np.ndarray, eps: float = 1e-5) -> np.ndarray:
        """Numerical gradient (central differences); used for Lipschitz checks."""
        n = self.n_atoms
        g = np.zeros(I.shape)
        for j in range(n):
            ej = np.zeros(n); ej[j] = eps
            g[..., j] = (self.energy(I + ej) - self.energy(I - ej)) / (2 * eps)
        return g


def make_random_horn_kb(n_atoms: int, n_clauses: int,
                        max_body: int = 3,
                        weight_scale: float = 1.0) -> HornKB:
    clauses = []
    for _ in range(n_clauses):
        body_cap = min(max_body, max(0, n_atoms - 1))
        body_len = int(RNG.integers(0, body_cap + 1))
        body = tuple(RNG.choice(n_atoms, size=body_len, replace=False).tolist())
        candidates = [j for j in range(n_atoms) if j not in body]
        head = int(RNG.choice(candidates)) if candidates else 0
        w = float(weight_scale * RNG.uniform(0.5, 1.5))
        clauses.append((head, body, w))
    return HornKB(n_atoms=n_atoms, clauses=clauses)


# ---------------------------------------------------------------------------
# Experiment 1: partition function via Monte Carlo & log Z bounds
# ---------------------------------------------------------------------------

def estimate_log_Z(kb: HornKB, n_samples: int = 200_000) -> tuple[float, float]:
    """Importance-sampling estimate of log Z = log integral exp(E) dI.

    Reference measure is uniform on [0,1]^n, so
        Z = E_uniform[exp(E(I))]
        log Z = log mean(exp(E(I))) computed via logsumexp for stability.
    Returns (log_Z_hat, standard_error).
    """
    I = RNG.uniform(0, 1, size=(n_samples, kb.n_atoms))
    e = kb.energy(I)
    log_Z = logsumexp(e) - np.log(n_samples)
    # Delta-method standard error on log Z
    w = np.exp(e - e.max())
    ess = (w.sum() ** 2) / (w ** 2).sum()
    se = 1.0 / np.sqrt(ess)
    return float(log_Z), float(se)


def theoretical_log_Z_bounds(kb: HornKB) -> tuple[float, float]:
    """Theoretical lower / upper bounds on log Z from Theorem 4.4 (volumetric)
    and Proposition 3.1 (basic): -W <= log Z <= W."""
    W = kb.W
    return -W, W


# ---------------------------------------------------------------------------
# Experiment 2: reflected Langevin diffusion & empirical mixing
# ---------------------------------------------------------------------------

def reflected_langevin(kb: HornKB,
                       n_chains: int = 200,
                       n_steps: int = 8000,
                       step_size: float = 5e-3,
                       I0: np.ndarray | None = None) -> np.ndarray:
    """Discretise the reflected overdamped Langevin SDE
        dI = grad E dt + sqrt(2) dB + dl,
    with reflection at the boundary of [0,1]^n.

    Returns the trajectory of shape (n_steps + 1, n_chains, n_atoms).
    """
    n = kb.n_atoms
    if I0 is None:
        I0 = RNG.uniform(0, 1, size=(n_chains, n))
    traj = np.empty((n_steps + 1, n_chains, n))
    traj[0] = I0
    sqrt2h = np.sqrt(2.0 * step_size)
    for t in range(n_steps):
        I = traj[t]
        g = kb.grad_energy(I)
        I_next = I + step_size * g + sqrt2h * RNG.standard_normal(size=I.shape)
        # Specular reflection at the cube boundary
        I_next = np.where(I_next < 0, -I_next, I_next)
        I_next = np.where(I_next > 1, 2.0 - I_next, I_next)
        # Clip residual escapes (pathological at corners)
        traj[t + 1] = np.clip(I_next, 0.0, 1.0)
    return traj


def empirical_mixing_time(traj: np.ndarray, tol: float = 0.05) -> int:
    """First step at which the empirical mean across chains stabilises
    within +/- tol of the long-run mean. Crude but matches the eta-norm
    convergence used in Cor. 7.5 of the manuscript."""
    mean_traj = traj.mean(axis=1)              # (n_steps + 1, n_atoms)
    long_run = mean_traj[-traj.shape[0] // 4:].mean(axis=0)
    deviation = np.linalg.norm(mean_traj - long_run, axis=1)
    converged = np.where(deviation < tol)[0]
    return int(converged[0]) if converged.size else traj.shape[0] - 1


def theoretical_mixing_time(kb: HornKB, eps: float = 0.05) -> float:
    """Cor. 7.5 of the manuscript:
        t_mix(eps) <= (e^{2W} / pi^2) * log(1 / eps)."""
    W = kb.W
    return float(np.exp(2 * W) * np.log(1.0 / eps) / np.pi ** 2)


# ---------------------------------------------------------------------------
# Experiment 3: Laplace approximation in the large-beta regime
# ---------------------------------------------------------------------------

def laplace_log_Z(kb: HornKB, beta: float, n_mc: int = 100_000) -> dict:
    """Test the leading-order Laplace asymptotic of Theorem 4.5.

    For a knowledge base with a unique vertex maximiser I^* and constant gradient
    g on the cell containing I^*, Theorem 4.5 predicts
        log Z(beta * w_hat) = beta * E_max - sum_j log(beta * |g_j|) + O(1)
    as beta -> infinity. We verify that the empirical log Z grows linearly
    in beta with slope E_max (the leading-order prediction).
    """
    n = kb.n_atoms
    w_norm = float(sum(abs(w) for _, _, w in kb.clauses))

    scaled_kb = HornKB(
        n_atoms=n,
        clauses=[(h, body, beta * w / w_norm) for h, body, w in kb.clauses],
    )
    log_Z_emp, _ = estimate_log_Z(scaled_kb, n_samples=n_mc)

    # Locate the vertex maximiser by exhaustive search over {0,1}^n
    if n <= 12:
        verts = np.array(np.unravel_index(np.arange(2 ** n),
                                          (2,) * n)).T.astype(float)
    else:
        verts = (RNG.uniform(0, 1, size=(2 ** 12, n)) > 0.5).astype(float)
    # Compute the energy at vertices for the *normalised* weight vector w_hat
    normalised_kb = HornKB(
        n_atoms=n,
        clauses=[(h, body, w / w_norm) for h, body, w in kb.clauses],
    )
    e_max_normalised = float(normalised_kb.energy(verts).max())

    # Leading-order Laplace prediction (Theorem 4.5):
    # log Z(beta * w_hat) ~ beta * E_max(w_hat) - n log(beta) + C
    log_Z_laplace_leading = beta * e_max_normalised
    return {
        "beta": beta,
        "E_max_normalised": e_max_normalised,
        "log_Z_empirical": log_Z_emp,
        "log_Z_leading_order": log_Z_laplace_leading,
        "log_Z_minus_leading": log_Z_emp - log_Z_laplace_leading,
        "expected_subleading_n_log_beta": -n * np.log(beta) if beta > 1 else 0.0,
    }


# ---------------------------------------------------------------------------
# Experiment 4: dimension scaling of the weight norm W
# ---------------------------------------------------------------------------

def scaling_study(n_grid: list[int],
                  density: float = 0.5,
                  n_repeats: int = 20) -> list[dict]:
    """For each n in n_grid, sample n_repeats Horn-clause KBs whose number
    of clauses scales as ceil(density * n * log n) and record the empirical
    W = ||w||_1. Connects to the dimension-scaling discussion of Section 9.4.
    """
    out = []
    for n in n_grid:
        m = max(1, int(np.ceil(density * n * np.log(max(n, 2)))))
        Ws = []
        for _ in range(n_repeats):
            kb = make_random_horn_kb(n_atoms=n, n_clauses=m, weight_scale=1.0)
            Ws.append(kb.W)
        out.append({
            "n": n,
            "m": m,
            "W_mean": float(np.mean(Ws)),
            "W_std": float(np.std(Ws)),
            "W_median": float(np.median(Ws)),
        })
    return out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run_all(output_dir: str = "results") -> dict:
    os.makedirs(output_dir, exist_ok=True)

    kb = make_random_horn_kb(n_atoms=6, n_clauses=10, weight_scale=0.5)
    W = kb.W

    # 1. Partition function bounds
    log_Z_emp, log_Z_se = estimate_log_Z(kb, n_samples=200_000)
    lo, hi = theoretical_log_Z_bounds(kb)
    partition = {
        "W": W, "log_Z_empirical": log_Z_emp, "log_Z_se": log_Z_se,
        "log_Z_lower_bound": lo, "log_Z_upper_bound": hi,
        "bound_satisfied": bool(lo <= log_Z_emp <= hi),
    }

    # 2. Mixing time
    traj = reflected_langevin(kb, n_chains=120, n_steps=4000, step_size=5e-3)
    t_emp = empirical_mixing_time(traj, tol=0.05)
    t_thm = theoretical_mixing_time(kb, eps=0.05)
    mixing = {
        "W": W,
        "empirical_mixing_steps": t_emp,
        "step_size": 5e-3,
        "empirical_mixing_time_continuous": t_emp * 5e-3,
        "theoretical_upper_bound": t_thm,
        "ratio_emp_over_thm": (t_emp * 5e-3) / max(t_thm, 1e-12),
    }

    # 3. Laplace asymptotic
    laplace = [laplace_log_Z(kb, beta=b, n_mc=30_000) for b in (1.0, 4.0, 16.0)]

    # 4. Scaling
    scaling = scaling_study(n_grid=[2, 4, 8, 16, 32, 64], density=0.5,
                            n_repeats=15)

    summary = {
        "partition": partition,
        "mixing": mixing,
        "laplace": laplace,
        "scaling": scaling,
    }
    with open(os.path.join(output_dir, "results.json"), "w") as f:
        json.dump(summary, f, indent=2)
    return summary


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir", default="results")
    args = p.parse_args()
    s = run_all(args.output_dir)
    print(json.dumps(s, indent=2))


if __name__ == "__main__":
    main()
