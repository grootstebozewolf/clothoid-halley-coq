using System.Runtime.CompilerServices;

namespace Clothoid.Halley;

/// <summary>
/// Result of a single solve: the arc length L (metres) and the number of
/// iterations consumed. Iterations == <see cref="ClothoidSolver.MaxIterDefault"/>
/// indicates non-convergence within the iteration budget.
/// </summary>
public readonly record struct SolverResult(double L, int Iterations);

/// <summary>
/// Halley's-method / Newton's-method solver for the clothoid
/// chord-length residual
///
///     f(L) = L^2 (P(L)^2 + Q(L)^2) - d^2,
///
/// where (P, Q, R, T, S2c, S2s) are the six moment integrals from
/// 32-point Gauss-Legendre quadrature on [0, 1] and d = |P1 - P0|.
/// Bit-identical reference: <c>python/clothoid_halleyL_bench.py</c>.
/// </summary>
public static class ClothoidSolver
{
    /// <summary>Default convergence tolerance.</summary>
    public const double TolDefault = 1e-13;

    /// <summary>Default iteration budget.</summary>
    public const int MaxIterDefault = 50;

    /// <summary>
    /// Solve f(L) = 0 by Halley's method (cubic convergence) with the
    /// four safety heuristics described in the paper.
    /// </summary>
    public static SolverResult SolveHalleyL(
        double[] p0, double[] p1, double k0, double k1,
        double tol = TolDefault, int maxIter = MaxIterDefault)
    {
        double cx = p1[0] - p0[0];
        double cy = p1[1] - p0[1];
        double d2 = cx * cx + cy * cy;
        double d  = Math.Sqrt(d2);
        if (d == 0.0) return new SolverResult(0.0, 0);

        double L = d;
        for (int it = 1; it <= maxIter; it++)
        {
            Moments(L, k0, k1, out double P, out double Q, out double R, out double T,
                                 out double S2c, out double S2s);
            double r2   = P * P + Q * Q;
            double qrpt = Q * R - P * T;
            double f    = L * L * r2 - d2;
            double fp   = 2.0 * L * r2 + 2.0 * L * L * qrpt;
            double fpp  = 2.0 * r2 + 8.0 * L * qrpt
                          + 2.0 * L * L * (R * R + T * T - P * S2c - Q * S2s);
            if (Math.Abs(f) < tol * Math.Max(d2, 1.0))
                return new SolverResult(L, it);

            double denom = 2.0 * fp * fp - f * fpp;
            if (Math.Abs(denom) < 1e-20 || fp <= 0.0)
            {
                L *= 1.5;
                continue;
            }
            double step  = 2.0 * f * fp / denom;
            double lNew  = L - step;
            if (lNew <= 0.0) lNew = 0.5 * L;
            L = lNew;
        }
        return new SolverResult(L, maxIter);
    }

    /// <summary>
    /// Solve f(L) = 0 by Newton's method (quadratic convergence).
    /// Same safety heuristics minus the Halley denominator guard.
    /// </summary>
    public static SolverResult SolveNewtonL(
        double[] p0, double[] p1, double k0, double k1,
        double tol = TolDefault, int maxIter = MaxIterDefault)
    {
        double cx = p1[0] - p0[0];
        double cy = p1[1] - p0[1];
        double d2 = cx * cx + cy * cy;
        double d  = Math.Sqrt(d2);
        if (d == 0.0) return new SolverResult(0.0, 0);

        double L = d;
        for (int it = 1; it <= maxIter; it++)
        {
            Moments(L, k0, k1, out double P, out double Q, out double R, out double T,
                                 out _, out _);
            double r2 = P * P + Q * Q;
            double f  = L * L * r2 - d2;
            double fp = 2.0 * L * r2 + 2.0 * L * L * (Q * R - P * T);
            if (Math.Abs(f) < tol * Math.Max(d2, 1.0))
                return new SolverResult(L, it);
            if (fp <= 0.0)
            {
                L *= 1.5;
                continue;
            }
            double lNew = L - f / fp;
            if (lNew < 0.5 * L) lNew = 0.5 * L;
            L = lNew;
        }
        return new SolverResult(L, maxIter);
    }

    /// <summary>
    /// Compute the six moment integrals
    ///   P, Q, R, T, S2c, S2s = ∫₀¹ {1, ψ, ψ²} · {cos, sin}(L·ψ(τ)) dτ
    /// with ψ(τ) = k0·τ + ½(k1-k0)·τ².
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Moments(double L, double k0, double k1,
        out double P, out double Q, out double R, out double T,
        out double S2c, out double S2s)
    {
        double half = 0.5 * (k1 - k0);
        double sP = 0, sQ = 0, sR = 0, sT = 0, sS2c = 0, sS2s = 0;
        for (int i = 0; i < GaussLegendre.N; i++)
        {
            double t   = GaussLegendre.T[i];
            double w   = GaussLegendre.W[i];
            double psi = k0 * t + half * t * t;
            double c   = Math.Cos(L * psi);
            double s   = Math.Sin(L * psi);
            sP   += w * c;
            sQ   += w * s;
            sR   += w * psi * c;
            sT   += w * psi * s;
            sS2c += w * psi * psi * c;
            sS2s += w * psi * psi * s;
        }
        P   = sP;
        Q   = sQ;
        R   = sR;
        T   = sT;
        S2c = sS2c;
        S2s = sS2s;
    }
}
