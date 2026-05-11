// SPDX-FileCopyrightText: 2026 Merkator Group
// SPDX-License-Identifier: LicenseRef-Merkator-Proprietary-NoAITraining

package com.grootstebozewolf.clothoid;

/**
 * Halley's / Newton's method for the clothoid chord-length residual
 *
 *     f(L) = L^2 (P(L)^2 + Q(L)^2) - d^2
 *
 * with six moment integrals (P, Q, R, T, S2c, S2s) from 32-point
 * Gauss-Legendre quadrature on [0, 1] and d = |P1 - P0|.
 *
 * Bit-identical reference: {@code python/clothoid_halleyL_bench.py}.
 */
public final class ClothoidSolver {

    /** Default convergence tolerance. */
    public static final double TOL_DEFAULT = 1e-13;

    /** Default iteration budget. */
    public static final int MAX_ITER_DEFAULT = 50;

    public record Result(double L, int iterations) {}

    private ClothoidSolver() {}

    /** Halley solve with default tolerance / iteration budget. */
    public static Result solveHalleyL(double[] p0, double[] p1, double k0, double k1) {
        return solveHalleyL(p0, p1, k0, k1, TOL_DEFAULT, MAX_ITER_DEFAULT);
    }

    /** Halley solve with explicit convergence parameters. */
    public static Result solveHalleyL(double[] p0, double[] p1, double k0, double k1,
                                      double tol, int maxIter) {
        double cx = p1[0] - p0[0];
        double cy = p1[1] - p0[1];
        double d2 = cx * cx + cy * cy;
        double d  = Math.sqrt(d2);
        if (d == 0.0) return new Result(0.0, 0);

        double[] m = new double[6];
        double L = d;
        for (int it = 1; it <= maxIter; it++) {
            moments(L, k0, k1, m);
            double P = m[0], Q = m[1], R = m[2], T = m[3], S2c = m[4], S2s = m[5];
            double r2   = P * P + Q * Q;
            double qrpt = Q * R - P * T;
            double f    = L * L * r2 - d2;
            double fp   = 2.0 * L * r2 + 2.0 * L * L * qrpt;
            double fpp  = 2.0 * r2 + 8.0 * L * qrpt
                        + 2.0 * L * L * (R * R + T * T - P * S2c - Q * S2s);
            if (Math.abs(f) < tol * Math.max(d2, 1.0)) {
                return new Result(L, it);
            }
            double denom = 2.0 * fp * fp - f * fpp;
            if (Math.abs(denom) < 1e-20 || fp <= 0.0) {
                L *= 1.5;
                continue;
            }
            double step = 2.0 * f * fp / denom;
            double lNew = L - step;
            if (lNew <= 0.0) lNew = 0.5 * L;
            L = lNew;
        }
        return new Result(L, maxIter);
    }

    public static Result solveNewtonL(double[] p0, double[] p1, double k0, double k1) {
        return solveNewtonL(p0, p1, k0, k1, TOL_DEFAULT, MAX_ITER_DEFAULT);
    }

    public static Result solveNewtonL(double[] p0, double[] p1, double k0, double k1,
                                      double tol, int maxIter) {
        double cx = p1[0] - p0[0];
        double cy = p1[1] - p0[1];
        double d2 = cx * cx + cy * cy;
        double d  = Math.sqrt(d2);
        if (d == 0.0) return new Result(0.0, 0);

        double[] m = new double[6];
        double L = d;
        for (int it = 1; it <= maxIter; it++) {
            moments(L, k0, k1, m);
            double P = m[0], Q = m[1], R = m[2], T = m[3];
            double r2 = P * P + Q * Q;
            double f  = L * L * r2 - d2;
            double fp = 2.0 * L * r2 + 2.0 * L * L * (Q * R - P * T);
            if (Math.abs(f) < tol * Math.max(d2, 1.0)) {
                return new Result(L, it);
            }
            if (fp <= 0.0) {
                L *= 1.5;
                continue;
            }
            double lNew = L - f / fp;
            if (lNew < 0.5 * L) lNew = 0.5 * L;
            L = lNew;
        }
        return new Result(L, maxIter);
    }

    /**
     * Compute the six moment integrals into out[0..5] = P, Q, R, T, S2c, S2s.
     * Java's `Math.fma` would be applicable here but not used; the simple form
     * matches the reference bit-for-bit on the test set.
     */
    private static void moments(double L, double k0, double k1, double[] out) {
        double half = 0.5 * (k1 - k0);
        double sP = 0, sQ = 0, sR = 0, sT = 0, sS2c = 0, sS2s = 0;
        for (int i = 0; i < GaussLegendre.N; i++) {
            double t   = GaussLegendre.T[i];
            double w   = GaussLegendre.W[i];
            double psi = k0 * t + half * t * t;
            double c   = Math.cos(L * psi);
            double s   = Math.sin(L * psi);
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
}
