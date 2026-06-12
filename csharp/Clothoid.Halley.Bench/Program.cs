// SPDX-FileCopyrightText: 2026 Merkator Group
// SPDX-License-Identifier: EUPL-1.2

using System.Diagnostics;
using System.IO;
using System.Text.Json;
using System.Text.Json.Serialization;
using Clothoid.Halley;

namespace Clothoid.Halley.Bench;

internal static class Program
{
    /// <summary>Warmup repetitions of the full corpus.</summary>
    private const int Warmups = 5;

    /// <summary>Measurement repetitions of the full corpus.</summary>
    private const int Repeats = 50;

    private static int Main()
    {
        var cases = LoadGoldenCases();
        var P0 = new double[cases.Length][];
        var P1 = new double[cases.Length][];
        var K0 = new double[cases.Length];
        var K1 = new double[cases.Length];
        for (int i = 0; i < cases.Length; i++)
        {
            P0[i] = cases[i].P0;
            P1[i] = cases[i].P1;
            K0[i] = cases[i].K0;
            K1[i] = cases[i].K1;
        }

        // Warmup (JIT tiering + cache fill)
        for (int w = 0; w < Warmups; w++)
        {
            _ = RunHalley(P0, P1, K0, K1);
            _ = RunNewton(P0, P1, K0, K1);
        }

        var halleyMs = new double[Repeats];
        var newtonMs = new double[Repeats];
        long halleyIters = 0, newtonIters = 0;
        for (int r = 0; r < Repeats; r++)
        {
            var sw = Stopwatch.StartNew();
            halleyIters = RunHalley(P0, P1, K0, K1);
            sw.Stop();
            halleyMs[r] = sw.Elapsed.TotalMilliseconds;

            sw.Restart();
            newtonIters = RunNewton(P0, P1, K0, K1);
            sw.Stop();
            newtonMs[r] = sw.Elapsed.TotalMilliseconds;
        }

        Array.Sort(halleyMs);
        Array.Sort(newtonMs);
        double halleyMedianMs = halleyMs[Repeats / 2];
        double newtonMedianMs = newtonMs[Repeats / 2];
        double halleyMinMs    = halleyMs[0];
        double newtonMinMs    = newtonMs[0];

        var result = new BenchResult(
            language:   "C# (.NET 8)",
            runtime:    Environment.Version.ToString(),
            cases:      cases.Length,
            halley_us:  1000.0 * halleyMedianMs / cases.Length,
            newton_us:  1000.0 * newtonMedianMs / cases.Length,
            halley_us_min: 1000.0 * halleyMinMs / cases.Length,
            newton_us_min: 1000.0 * newtonMinMs / cases.Length,
            halley_iter_mean: (double)halleyIters / cases.Length,
            newton_iter_mean: (double)newtonIters / cases.Length
        );
        Console.WriteLine(JsonSerializer.Serialize(result, new JsonSerializerOptions { WriteIndented = true }));
        return 0;
    }

    private static long RunHalley(double[][] P0, double[][] P1, double[] K0, double[] K1)
    {
        long iters = 0;
        for (int i = 0; i < P0.Length; i++)
        {
            var r = ClothoidSolver.SolveHalleyL(P0[i], P1[i], K0[i], K1[i]);
            iters += r.Iterations;
        }
        return iters;
    }

    private static long RunNewton(double[][] P0, double[][] P1, double[] K0, double[] K1)
    {
        long iters = 0;
        for (int i = 0; i < P0.Length; i++)
        {
            var r = ClothoidSolver.SolveNewtonL(P0[i], P1[i], K0[i], K1[i]);
            iters += r.Iterations;
        }
        return iters;
    }

    private static Case[] LoadGoldenCases()
    {
        var path = FindRepoRelative("data/golden_vectors.json");
        using var fs = File.OpenRead(path);
        var doc = JsonSerializer.Deserialize<File_>(fs, JsonOpts)
            ?? throw new InvalidDataException("failed to parse golden_vectors.json");
        return doc.Cases.ToArray();
    }

    private static readonly JsonSerializerOptions JsonOpts = new()
    {
        PropertyNamingPolicy = JsonNamingPolicy.SnakeCaseLower,
    };

    private static string FindRepoRelative(string relative)
    {
        var dir = new DirectoryInfo(AppContext.BaseDirectory);
        while (dir != null)
        {
            var candidate = Path.Combine(dir.FullName, relative.Replace('/', Path.DirectorySeparatorChar));
            if (File.Exists(candidate)) return candidate;
            dir = dir.Parent;
        }
        throw new FileNotFoundException($"could not locate '{relative}'");
    }

    private sealed class File_
    {
        [JsonPropertyName("cases")] public List<Case> Cases { get; set; } = new();
    }

    private sealed class Case
    {
        [JsonPropertyName("P0")] public double[] P0 { get; set; } = System.Array.Empty<double>();
        [JsonPropertyName("P1")] public double[] P1 { get; set; } = System.Array.Empty<double>();
        [JsonPropertyName("k0")] public double K0 { get; set; }
        [JsonPropertyName("k1")] public double K1 { get; set; }
    }

    // ReSharper disable NotAccessedPositionalProperty.Global
    private sealed record BenchResult(
        string language,
        string runtime,
        int    cases,
        double halley_us,
        double newton_us,
        double halley_us_min,
        double newton_us_min,
        double halley_iter_mean,
        double newton_iter_mean);
}
