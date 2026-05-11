// SPDX-FileCopyrightText: 2026 Merkator Group
// SPDX-License-Identifier: LicenseRef-Merkator-Proprietary-NoAITraining

/**
 * Halley's / Newton's method for the clothoid chord-length residual
 *
 *     f(L) = L^2 (P(L)^2 + Q(L)^2) - d^2
 *
 * with six moment integrals (P, Q, R, T, S2c, S2s) from 32-point
 * Gauss-Legendre quadrature on [0, 1] and d = |P1 - P0|.
 *
 * Bit-identical reference: python/clothoid_halleyL_bench.py
 */
import { GL_N, GL_T, GL_W } from "./gaussLegendre.js";

export interface SolverResult {
    readonly L:          number;
    readonly iterations: number;
}

export const TOL_DEFAULT     = 1e-13;
export const MAX_ITER_DEFAULT = 50;

export function solveHalleyL(
    p0:       readonly [number, number],
    p1:       readonly [number, number],
    k0:       number,
    k1:       number,
    tol:      number = TOL_DEFAULT,
    maxIter:  number = MAX_ITER_DEFAULT,
): SolverResult {
    const cx = p1[0] - p0[0];
    const cy = p1[1] - p0[1];
    const d2 = cx * cx + cy * cy;
    const d  = Math.sqrt(d2);
    if (d === 0) return { L: 0, iterations: 0 };

    let L = d;
    const m = new Float64Array(6);
    for (let it = 1; it <= maxIter; it++) {
        moments(L, k0, k1, m);
        const P = m[0], Q = m[1], R = m[2], T = m[3], S2c = m[4], S2s = m[5];
        const r2   = P * P + Q * Q;
        const qrpt = Q * R - P * T;
        const f    = L * L * r2 - d2;
        const fp   = 2 * L * r2 + 2 * L * L * qrpt;
        const fpp  = 2 * r2 + 8 * L * qrpt
                   + 2 * L * L * (R * R + T * T - P * S2c - Q * S2s);
        if (Math.abs(f) < tol * Math.max(d2, 1)) {
            return { L, iterations: it };
        }
        const denom = 2 * fp * fp - f * fpp;
        if (Math.abs(denom) < 1e-20 || fp <= 0) {
            L *= 1.5;
            continue;
        }
        const step = (2 * f * fp) / denom;
        let lNew = L - step;
        if (lNew <= 0) lNew = 0.5 * L;
        L = lNew;
    }
    return { L, iterations: maxIter };
}

export function solveNewtonL(
    p0:       readonly [number, number],
    p1:       readonly [number, number],
    k0:       number,
    k1:       number,
    tol:      number = TOL_DEFAULT,
    maxIter:  number = MAX_ITER_DEFAULT,
): SolverResult {
    const cx = p1[0] - p0[0];
    const cy = p1[1] - p0[1];
    const d2 = cx * cx + cy * cy;
    const d  = Math.sqrt(d2);
    if (d === 0) return { L: 0, iterations: 0 };

    let L = d;
    const m = new Float64Array(6);
    for (let it = 1; it <= maxIter; it++) {
        moments(L, k0, k1, m);
        const P = m[0], Q = m[1], R = m[2], T = m[3];
        const r2 = P * P + Q * Q;
        const f  = L * L * r2 - d2;
        const fp = 2 * L * r2 + 2 * L * L * (Q * R - P * T);
        if (Math.abs(f) < tol * Math.max(d2, 1)) {
            return { L, iterations: it };
        }
        if (fp <= 0) {
            L *= 1.5;
            continue;
        }
        let lNew = L - f / fp;
        if (lNew < 0.5 * L) lNew = 0.5 * L;
        L = lNew;
    }
    return { L, iterations: maxIter };
}

function moments(L: number, k0: number, k1: number, out: Float64Array): void {
    const half = 0.5 * (k1 - k0);
    let sP = 0, sQ = 0, sR = 0, sT = 0, sS2c = 0, sS2s = 0;
    for (let i = 0; i < GL_N; i++) {
        const t   = GL_T[i]!;
        const w   = GL_W[i]!;
        const psi = k0 * t + half * t * t;
        const c   = Math.cos(L * psi);
        const s   = Math.sin(L * psi);
        sP   += w * c;
        sQ   += w * s;
        sR   += w * psi * c;
        sT   += w * psi * s;
        sS2c += w * psi * psi * c;
        sS2s += w * psi * psi * s;
    }
    out[0] = sP;
    out[1] = sQ;
    out[2] = sR;
    out[3] = sT;
    out[4] = sS2c;
    out[5] = sS2s;
}
